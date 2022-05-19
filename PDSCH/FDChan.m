classdef FDChan < matlab.System
    % Frequency-domain multipath channel
    properties
        % Configuration
        carrierConfig;   % Carrier configuration
        
        % Path parameters
        gain;  % Path gain in dB
        dly;   % Delay of each path in seconds
        aoaAz, aoaEl; % Angle of arrival of each path in degrees      
        fd;    % Doppler shift for each path
        
        rxVel = [30,0,0]';  % Mobile velocity vector in m/s
        fc = 28e9;    % Carrier freq in Hz
        
        gainComplex;  % Complex gain of each path
        
        % SNR parameters
        Etx = 1;       % average energy per PDSCH symbol 
        EsN0Avg;  % Avg SNR per RX symbol in dB
       
        % Symbol times
        symStart;  % symStart(i) = start of symbol i relative to subframe

        
        waveformConfig;
                     
    end
    methods
        function obj = FDChan(carrierConfig, varargin)
            % Constructor
            
          
            obj.carrierConfig = carrierConfig;

                                 
           
            if nargin >= 1
                obj.set(varargin{:});
            end
            
          
            
            
            obj.gainComplex = db2pow(obj.gain).*exp(-1i*2*pi*rand([2 1])) ;          
           
          

            %assuming only one frame of reference
            
            % Computing  unit vector in direction of each path

            [x,y,z] = sph2cart(deg2rad(obj.aoaAz),deg2rad(obj.aoaEl),1);
            u = [x;y;z];

            obj.fd = obj.rxVel(1,1).*cos(obj.aoaAz)*obj.fc/physconst('LightSpeed');
            
    
            
            obj.symStart = zeros(1,14*8);
            total_sym = 0;
            a = nrOFDMInfo(carrierConfig);

            for i = 1:112
            
            total_sym = total_sym + a.SymbolLengths(i);
            obj.symStart(i) = total_sym/a.SampleRate;

            end
                                                   
        end
        
        
    end
    methods (Access = protected)
        
        
        function [rxGrid, chanGrid, noiseVar] = stepImpl(obj, txGrid, sfNum, slotNum)
            % Applies a frequency domain channel and noise
           
           

            t_frame = obj.symStart + sfNum*obj.symStart(end );
            t_slot = t_frame(slotNum*14 +1 :(slotNum+1)*14);
            gainLin = (sum(sqrt(abs(obj.gainComplex))))^2;
            noiseVar = obj.Etx*gainLin/db2pow(obj.EsN0Avg);
            
            chanGrid = zeros(792,14);

           for n = 1 : 792
                for j = 1 :14
            chanGrid(n,j) = sum(obj.gainComplex.*exp(2*pi*1i*(t_slot(j)*obj.fd + n*120e3*obj.dly)),1);
                end
           end
           
            rxGrid = chanGrid.*txGrid + sqrt(noiseVar).*rand();
   
            
                
           
            
 
        end
        
    end
end