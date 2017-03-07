function acqResults = acquisition(longSignal, settings)

%longSignal:�����Ҳ����ݣ���ӳ��ͼ����Ϊ���Ҳ�����������
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 11 ms of raw signal from the front-end 
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the 
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number. 
 
%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis and Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on Peter Rinder and Nicolaj Bertelsen
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $

%% Initialization =========================================================

% Find number of samples per spreading code
%ÿ����Ƭ�Ĳ�����
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

% Create two 1msec vectors of data to correlate with and one with zero DC
%��������1ms���ź�
signal1 = longSignal(1 : samplesPerCode);
signal2 = longSignal(samplesPerCode+1 : 2*samplesPerCode);

%��ֱ���ź�
signal0DC = longSignal - mean(longSignal); 

% Find sampling period
%��������
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave 
%ȥ��Ƶ��f����λ��
phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band (500Hz steps)
%Ƶ����������500Hz����
numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;

% Generate all C/A codes and sample them according to the sampling freq.
%�������C/A���
caCodesTable = makeCaTable(settings);


%--- Initialize arrays to speed up the code -------------------------------
% Search results of all frequency bins and code shifts (for one satellite)
%results��ĳ��������ĳ���ض�Ƶ���������¶�Ӧ����ؽ������һ��Ƶ������������*������������еľ���
results     = zeros(numberOfFrqBins, samplesPerCode);

% Carrier frequencies of the frequency bins
frqBins     = zeros(1, numberOfFrqBins);


%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 37);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 37);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 37);

fprintf('(');

% Perform search for all listed PRN numbers ...
%�������е����Ǻ�
for PRN = settings.acqSatelliteList
%�ο��飺��������GPS��٤���Խ��ջ���74ҳͼ6.8
%% Correlate signals ======================================================   
    %--- Perform DFT of C/A code ------------------------------------------
    %C/A�븵��Ҷ�任����ȡ����
    caCodeFreqDom = conj(fft(caCodesTable(PRN, :)));
    %--- Make the correlation for whole frequency band (for all freq. bins)
    %������ǰ�������е�Ƶ�㣬500Hz����
    for frqBinIndex = 1:numberOfFrqBins
    
        %--- Generate carrier wave frequency grid (0.5kHz step) -----------
        %����ÿ��Ƶ���Ƶ��ֵ
        frqBins(frqBinIndex) = settings.IF - ...
                               (settings.acqSearchBand/2) * 1000 + ...
                               0.5e3 * (frqBinIndex - 1);

        %--- Generate local sine and cosine -------------------------------
        %���������ز��źţ�I·��Q·
        sinCarr = sin(frqBins(frqBinIndex) * phasePoints);
        cosCarr = cos(frqBins(frqBinIndex) * phasePoints);

        %--- "Remove carrier" from the signal -----------------------------
        %ȥ�ز���������λ����
        I1      = sinCarr .* signal1;
        Q1      = cosCarr .* signal1;
        I2      = sinCarr .* signal2;
        Q2      = cosCarr .* signal2;

        %--- Convert the baseband signal to frequency domain --------------
        %����Ҷ�任��׼����ȡ�����C/A�����
        IQfreqDom1 = fft(I1 + j*Q1);
        IQfreqDom2 = fft(I2 + j*Q2);

        %--- Multiplication in the frequency domain (correlation in time
        %domain)
        %����Ҷ�任���ȥ�ز��ź���ȡ������C/A�����
        convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
        convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;

        %--- Perform inverse DFT and store correlation results ------------
        %����˽������Ҷ���任��ȡƽ��
        acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;
        acqRes2 = abs(ifft(convCodeIQ2)) .^ 2;
        
        %--- Check which msec had the greater power and save that, will
        %"blend" 1st and 2nd msec but will correct data bit issues
        %�ҳ���ǰ2ms�ź������ֵ�ϴ���źţ��洢��Ƶ���Ӧ�����ֵ
        if (max(acqRes1) > max(acqRes2))
            results(frqBinIndex, :) = acqRes1;
        else
            results(frqBinIndex, :) = acqRes2;
        end
        
    end % frqBinIndex = 1:numberOfFrqBins

%% Look for correlation peaks in the results ==============================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
    
    %--- Find the correlation peak and the carrier frequency --------------
    %max(results, [], 2)��ʾresults��ÿ�е�����ֵ
    %�ⲿ��max��ʾresult�����ֵ����¼���ֵ�Լ�����λ�ã���Ƶ�ʵ�������
    %peakSize:���ֵ
    %frequencyBinIndex��Ƶ�ʵ�����
    [peakSize frequencyBinIndex] = max(max(results, [], 2));
    %--- Find code phase of the same correlation peak ---------------------
    %�ҳ����ֵ���ڵ�����λλ��
    [peakSize codePhase] = max(max(results));
    
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    %ÿ����Ƭ�Ĳ�����
    samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
    %�ҵ���ֵ���ڵ�����λ����ǰ�����ȥһ����Ƭ�ľ��룬�������Χ��Ѱ��ͬһƵ�ʵĵڶ����ֵ
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;

    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    %�ڷ�ֵ����λ����ǰ��һ����Ƭ��Ѱ�ҵڶ����ֵ
    
    %�˴���Դ����ΪС��2�����ڱ��˸�ΪС��1.����Ϊ��С��2��������Ϊ�Ǵ�������
    %�����·���ź����ӳٰ˷�֮�˸���Ƭʱ��codePhase = 9ʱ�������ˡ�Index exceeds matrix dimensions������
    %if excludeRangeIndex1 < 2
    if excludeRangeIndex1 < 1
        codePhaseRange = excludeRangeIndex2 : ...
                         (samplesPerCode + excludeRangeIndex1);
                         
    elseif excludeRangeIndex2 >= samplesPerCode
        codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                         excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
                          excludeRangeIndex2 : samplesPerCode];
    end

    %--- Find the second highest correlation peak in the same freq. bin ---
    %secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));
    [secondPeakSize secondCodePhase] = max(results(frequencyBinIndex, codePhaseRange));
    
    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(PRN) = peakSize/secondPeakSize;
    if frequencyBinIndex == 15 && PRN == 2
        figure;
        plot(results(frequencyBinIndex, :));
        fprintf('results1~20');
        results(frequencyBinIndex, 1:20)
        fprintf('results16348~16368');
        results(frequencyBinIndex, samplesPerCode - 19: samplesPerCode)
        codePhase
        secondCodePhase
        peakSize
        secondPeakSize
        fprintf('peakSize/secondPeakSize');
        acqResults.peakMetric(PRN)
    end
    % If the result is above threshold, then there is a signal ...
    %������ֵ��ڶ���ֵ��ֵ������ֵ������������
    if (peakSize/secondPeakSize) > settings.acqThreshold

%% Fine resolution frequency search =======================================
%ϸ������Ƶ��
        %--- Indicate PRN number of the detected signal -------------------
        fprintf('%02d ', PRN);
        
        %--- Generate 10msec long C/A codes sequence for given PRN --------
        caCode = generateCAcode(PRN);
        
        codeValueIndex = floor((ts * (1:10*samplesPerCode)) / ...
                               (1/settings.codeFreqBasis));
                           
        longCaCode = caCode((rem(codeValueIndex, 1023) + 1));
    
        %--- Remove C/A code modulation from the original signal ----------
        % (Using detected C/A code phase)
        xCarrier = ...
            signal0DC(codePhase:(codePhase + 10*samplesPerCode-1)) ...
            .* longCaCode;
        
        %--- Find the next highest power of two and increase by 8x --------
        fftNumPts = 8*(2^(nextpow2(length(xCarrier))));
        
        %--- Compute the magnitude of the FFT, find maximum and the
        %associated carrier frequency 
        fftxc = abs(fft(xCarrier, fftNumPts)); 
        
        uniqFftPts = ceil((fftNumPts + 1) / 2);
        [fftMax, fftMaxIndex] = max(fftxc(5 : uniqFftPts-5));
        
        fftFreqBins = (0 : uniqFftPts-1) * settings.samplingFreq/fftNumPts;
        
        %--- Save properties of the detected satellite signal -------------
        acqResults.carrFreq(PRN)  = fftFreqBins(fftMaxIndex);
        acqResults.codePhase(PRN) = codePhase;
        
        %��������ˣ��ͻ�ͼ
        fprintf('�����ˣ�');
        codePhase
        secondCodePhase
        frequencyBinIndex
        figure;
        %plot(1:samplesPerCode,results(frequencyBinIndex, :));
        plot(1:30,results(frequencyBinIndex, 1:30));
        fprintf('results1~5;16364~16368');
        results(frequencyBinIndex, 1:5)
        results(frequencyBinIndex, samplesPerCode - 4 : samplesPerCode)
        %������ͼ
        figure
        bar(1:30,results(frequencyBinIndex, 1:30));
        
        
        frequencyBinIndex
        samplesPerCodeChip
    
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
        
    
end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
