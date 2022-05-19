classdef NRgNBTxFD < matlab.System
    % 5G NR gNB transmitter class implemented in frequency domain
    properties
        % Configuration
        carrierConfig;   % Carrier configuration
        pdschConfig;     % PDSCH configuration
                                               
        % Transport block data for last transmission
        targetCodeRate = 490/1024;  % Target code rate
        trBlkSizes;                 % Transport block size
        
        % Transmitted data in last slots
        txBits;         % TX bits
               
        % DLSCH encoder
        encDLSCH;        
        
        
    end
    methods
        function obj = NRgNBTxFD(carrierConfig, pdschConfig, ...
                varargin)
            % Constructor
            
            % Saving  the carrier and PDSCH configuration
            obj.carrierConfig = carrierConfig;
            obj.pdschConfig = pdschConfig;
                                                           
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
            
            % Creating  DLSCH encoder system object
            obj.encDLSCH = nrDLSCH('MultipleHARQProcesses', false, ...
                'TargetCodeRate', obj.targetCodeRate);   
                             
        end
    end
    methods (Access = protected)
               
        function txGrid = stepImpl(obj)
            % step implementation. Creates one slot of samples for each
            % component carrier
            
            
        
            txGrid = nrResourceGrid(obj.carrierConfig, ...
                obj.pdschConfig.NumLayers);
            
                       
            % Getting  indices on where the PDSCH is allocated
            [pdschInd,pdschInfo] = nrPDSCHIndices(obj.carrierConfig, obj.pdschConfig);
            
            % Computing  the extra overhead from the PT-RS
            Xoh_PDSCH = 6*obj.pdschConfig.EnablePTRS;     
            
            % Calculating  the transport block size based on the PDSCH
            % allocation and target code rate
            obj.trBlkSizes = nrTBS(obj.pdschConfig.Modulation,obj.pdschConfig.NumLayers,...
                numel(obj.pdschConfig.PRBSet),pdschInfo.NREPerPRB,...
                obj.targetCodeRate,Xoh_PDSCH);
            
            % Generating  random bits for each codeword and setting  the transport
            % block
            obj.txBits = cell(obj.pdschConfig.NumCodewords, 1);
            for icw = 1:obj.pdschConfig.NumCodewords

                % Creating random bits
                obj.txBits{icw} = randi([0 1], obj.trBlkSizes(icw), 1);
                
                % Encoding  data                
                obj.encDLSCH.setTransportBlock(obj.txBits{icw},icw-1);
            end
            
            % Encoding  the DL-SCH transport blocks.  
            rv = 0;
            codedTrBlock = obj.encDLSCH(obj.pdschConfig.Modulation, ...
                obj.pdschConfig.NumLayers, pdschInfo.G, rv);
            
            % Modulateing the PDSCH modulation
            pdschSymbols = nrPDSCH(obj.carrierConfig, obj.pdschConfig, ...
                codedTrBlock);
          
            % Mapping  the modulated symbols to the OFDM grid
            txGrid(pdschInd) = pdschSymbols;            
                                    
        end
        
    end
end