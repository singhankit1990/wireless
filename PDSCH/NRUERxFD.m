classdef NRUERxFD < matlab.System
    % 5G NR UR receiver class implemented in frequency domain
    properties
        % Configuration
        carrierConfig;   % Carrier configuration
        pdschConfig;     % Default PDSCH config
        waveformConfig;  % Waveform config
        
        % OFDM grid
        rxGrid;
        
        % Transport block data for last transmission
        targetCodeRate = 490/1024;  % Target code rate
        trBlkSizes;                 % Transport block size
        
        % Received data in last slots
        pdschEq;       % Equalized PDSCH symbols
        rxBits;        % RX bits
        
        % DLSCH decoder
        decDLSCH;
        
        
    end
    methods
        function obj = NRUERxFD(carrierConfig, pdschConfig, ...
                varargin)
            % Constructor
            
            % Saving the carrier and PDSCH configuration
            obj.carrierConfig = carrierConfig;
            obj.pdschConfig = pdschConfig;
            
            % Creating the waveform configuration from the carrier
            % configuration
            obj.waveformConfig = nrOFDMInfo(obj.carrierConfig);
            
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
            
            % Creating DLSCH decoder
            obj.decDLSCH = nrDLSCHDecoder('MultipleHARQProcesses', false, ...
                'TargetCodeRate', obj.targetCodeRate, ...
                'LDPCDecodingAlgorithm', 'Layered belief propagation');
            
        end
    end
    methods (Access = protected)
        
        
        function stepImpl(obj, rxGrid, chanGrid, noiseVar)
           
            [pdschInd,pdschInfo] = nrPDSCHIndices(...
                obj.carrierConfig, obj.pdschConfig);
            [pdschRx, pdschHest] = nrExtractResources(pdschInd, rxGrid,...
                chanGrid);
            
            % eqalizing           
            [obj.pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseVar);
            
                        
            % getting the LLRs 
            [dlschLLRs,rxSym] = nrPDSCHDecode(obj.carrierConfig,obj.pdschConfig,obj.pdschEq,noiseVar);   
            
          
            csi = nrLayerDemap(csi); % CSI layer demapping
            numCW = length(csi);
            for cwIdx = 1:numCW
                Qm = length(dlschLLRs{cwIdx})/length(rxSym{cwIdx}); % bits per symbol
                csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);   % expand by each bit per symbol
                dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   % scale
            end
                    
            % Computing the extra overhead from the PT-RS
            Xoh_PDSCH = 6*obj.pdschConfig.EnablePTRS;
            
            % Calculating the transport block size based on the PDSCH
            % allocation and target code rate
            obj.trBlkSizes = nrTBS(obj.pdschConfig.Modulation,obj.pdschConfig.NumLayers,...
                numel(obj.pdschConfig.PRBSet),pdschInfo.NREPerPRB,...
                obj.targetCodeRate,Xoh_PDSCH);
            obj.decDLSCH.TransportBlockLength = obj.trBlkSizes;
                         
            % Reset the soft buffer
            harqId = 0;
            obj.decDLSCH.resetSoftBuffer(harqId);
            
      
            rv = 0;
            obj.rxBits = obj.decDLSCH(dlschLLRs,'16QAM',1,rv);
            
        end
        
    end
end