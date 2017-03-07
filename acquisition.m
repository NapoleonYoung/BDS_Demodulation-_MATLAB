function acqResults = acquisition(longSignal, settings)

%longSignal:是正弦波数据，反映在图表上为正弦波采样点连线
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
%每个码片的采样数
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

% Create two 1msec vectors of data to correlate with and one with zero DC
%创建两个1ms的信号
signal1 = longSignal(1 : samplesPerCode);
signal2 = longSignal(samplesPerCode+1 : 2*samplesPerCode);

%零直流信号
signal0DC = longSignal - mean(longSignal); 

% Find sampling period
%采样周期
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave 
%去掉频率f的相位点
phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band (500Hz steps)
%频率搜索数，500Hz步进
numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;

% Generate all C/A codes and sample them according to the sampling freq.
%采样后的C/A码表
caCodesTable = makeCaTable(settings);


%--- Initialize arrays to speed up the code -------------------------------
% Search results of all frequency bins and code shifts (for one satellite)
%results：某颗卫星在某个特定频率搜索点下对应的相关结果，是一个频率搜索数个行*码采样点数个列的矩阵
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
%搜索所有的卫星号
for PRN = settings.acqSatelliteList
%参考书：软件定义的GPS和伽利略接收机第74页图6.8
%% Correlate signals ======================================================   
    %--- Perform DFT of C/A code ------------------------------------------
    %C/A码傅立叶变换，并取共轭
    caCodeFreqDom = conj(fft(caCodesTable(PRN, :)));
    %--- Make the correlation for whole frequency band (for all freq. bins)
    %搜索当前卫星所有的频点，500Hz步进
    for frqBinIndex = 1:numberOfFrqBins
    
        %--- Generate carrier wave frequency grid (0.5kHz step) -----------
        %产生每个频点的频率值
        frqBins(frqBinIndex) = settings.IF - ...
                               (settings.acqSearchBand/2) * 1000 + ...
                               0.5e3 * (frqBinIndex - 1);

        %--- Generate local sine and cosine -------------------------------
        %产生本地载波信号，I路与Q路
        sinCarr = sin(frqBins(frqBinIndex) * phasePoints);
        cosCarr = cos(frqBins(frqBinIndex) * phasePoints);

        %--- "Remove carrier" from the signal -----------------------------
        %去载波，不必相位对齐
        I1      = sinCarr .* signal1;
        Q1      = cosCarr .* signal1;
        I2      = sinCarr .* signal2;
        Q2      = cosCarr .* signal2;

        %--- Convert the baseband signal to frequency domain --------------
        %傅立叶变换，准备与取共轭的C/A码相乘
        IQfreqDom1 = fft(I1 + j*Q1);
        IQfreqDom2 = fft(I2 + j*Q2);

        %--- Multiplication in the frequency domain (correlation in time
        %domain)
        %傅立叶变换后的去载波信号与取共轭后的C/A码相乘
        convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
        convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;

        %--- Perform inverse DFT and store correlation results ------------
        %将相乘结果傅立叶反变换，取平方
        acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;
        acqRes2 = abs(ifft(convCodeIQ2)) .^ 2;
        
        %--- Check which msec had the greater power and save that, will
        %"blend" 1st and 2nd msec but will correct data bit issues
        %找出当前2ms信号中相关值较大的信号，存储该频点对应的相关值
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
    %max(results, [], 2)表示results中每行的最大的值
    %外部的max表示result中最大值，记录最大值以及所在位置（即频率点索引）
    %peakSize:最大值
    %frequencyBinIndex：频率点索引
    [peakSize frequencyBinIndex] = max(max(results, [], 2));
    %--- Find code phase of the same correlation peak ---------------------
    %找出最大值所在的码相位位置
    [peakSize codePhase] = max(max(results));
    
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    %每个码片的采样数
    samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
    %找到峰值所在的码相位处，前后各减去一个码片的距离，在这个范围内寻找同一频率的第二大峰值
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;

    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    %在峰值码相位索引前后一个码片外寻找第二大峰值
    
    %此处，源代码为小于2，现在本人改为小于1.至于为何小于2，本人认为是代码有误，
    %例如多路径信号是延迟八分之八个码片时，codePhase = 9时，出现了“Index exceeds matrix dimensions”错误
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
    %如果最大值与第二大值比值超出阈值，即捕获到卫星
    if (peakSize/secondPeakSize) > settings.acqThreshold

%% Fine resolution frequency search =======================================
%细化采样频率
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
        
        %如果捕获到了，就画图
        fprintf('捕获到了：');
        codePhase
        secondCodePhase
        frequencyBinIndex
        figure;
        %plot(1:samplesPerCode,results(frequencyBinIndex, :));
        plot(1:30,results(frequencyBinIndex, 1:30));
        fprintf('results1~5;16364~16368');
        results(frequencyBinIndex, 1:5)
        results(frequencyBinIndex, samplesPerCode - 4 : samplesPerCode)
        %画条形图
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
