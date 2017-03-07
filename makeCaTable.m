function [ sampledCATable ] = makeCaTable( settings )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%����һ�����ڵ�C/A���ܵĲ�������
samplesPerCode = round(settings.samplingFreq /(settings.codeFreqBasis / settings.codeLength));
%����Ĳ���C/A����
sampledCATable = zeros(37, samplesPerCode);

%����һ�ε�ʱ�䳤��
ts = 1 / settings.samplingFreq;

%һ����Ƭ��ʱ�䳤��
tc = 1 / settings.codeFreqBasis;

%PRN = 37ԭ��37������
for PRN = 1 : 37
    
    %����C/A��
    CACode = generateCAcode(PRN);
    
    %-----���ֻ�C/A��-----
    
    %��������C/A���λ������
    codeValueIndex = ceil((ts / tc) * (1 : samplesPerCode));
    
    %���һ�����������ɵ�λ��һ����C/A���һ��ֵ
    codeValueIndex(end) = settings.codeLength;
    
    sampledCATable(PRN, :) =  CACode(codeValueIndex);
end
    
end

