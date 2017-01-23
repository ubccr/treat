%% Copyright 2015 TREAT Authors. All rights reserved.
%%
%% This file is part of TREAT.
%%
%% TREAT is free software: you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation, either version 3 of the License, or
%% (at your option) any later version.
%%
%% TREAT is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with TREAT.  If not, see <http://www.gnu.org/licenses/>.
%%
%%-------------------------------------------------------------------------------------------
%%
%%----------------------------------------------------------------------%%
%   Statistical test for significant editing site selection.
%   Author: Runpu Chen. runpuche@buffalo.edu
%   Update date: Jan., 8, 2017
%%----------------------------------------------------------------------%%

clear all
close all
sourceFile = '/file_1013/';                           % the local director of the source files.
resultFile = '/results_1013/';                        % the local director of the result files.
allFiles = dir( strcat(pwd, sourceFile));
allNames_temp = { allFiles.name };
allNames = allNames_temp(4:end);
clear allFiles
num_file = length(allNames);
smallVal = 0.001;
for i_file = 1:length(allNames)
disp(strcat('Processing the number',{' '}, num2str(i_file),{' '}, 'file...'));
    filename = strcat(pwd, sourceFile, allNames{i_file});
    fid = fopen(filename);
    header = textscan(fid,'%s',1);
    header = strsplit(char(header{1}), ',');
    fclose(fid);
    header = header(1:end);
    data = csvread(filename, 1, 0);                  % read in data
    id = data(:, 1);
    data = data(:, 2:end);
    num_row = size(data, 1);
    num_column = size(data, 2);
    inducedNum = 2;                                 % the last two columns are induced
    
    inducedData = data(:, end-inducedNum+1:end);
    uninducedData = data(:, 1:end-inducedNum);
    headerUninduced = header(1:end-inducedNum);
    headerInduced = header(end-inducedNum+1:end);
    
    p = 0*inducedData; q = p;
    for irow = 1:num_row
        control =  uninducedData(irow, :);
        meanControl = mean(control);
        stdControl = std(control);
        for j = 1:inducedNum
            prob = normcdf(inducedData(irow,j),meanControl,stdControl);
            if prob >= 0.5
                p(irow, j) = 2*(1-prob);
            else p(irow,j) = 2*prob;
            end
            if meanControl == 0
                p(irow, j) = NaN;
            end
            if meanControl == 0 & inducedData(irow,j) > smallVal
                p(irow, j) = 0;
            end
            
        end
        for j = 1:inducedNum
            q(:, j) = mafdr(p(:, j),'BHFDR', 'true');
        end
    end
    
    headerFinal = headerUninduced;
    dataFinal = [id, uninducedData];
    for j = 1:inducedNum
        headerFinal = [headerFinal, headerInduced{j}, 'pVal', 'qVal'];
        dataFinal = [dataFinal, inducedData(:, j), p(:, j), q(:, j)];
        
    end
    resultFileName = strcat(pwd, resultFile,'sepResults/',allNames{i_file} , 'SepResult.txt');
    fid = fopen(resultFileName, 'w');           % write the result files.
    for j = 1:length(headerFinal)
        fprintf(fid, '%s\t', char(headerFinal{j})); %write header.
    end
    fprintf(fid, '\n');
    for i = 1:num_row
        fprintf(fid, '%d\t', dataFinal(i, 1));
        for iColumn = 2:length(headerFinal)
            fprintf(fid, '%8.10f\t', dataFinal(i, iColumn)); %write data.
        end
        fprintf(fid, '\n');
        
    end
    fclose(fid);
    %% for combined analyze
    p = zeros(num_row, 1); q = p;
    for irow = 1:num_row
        control =  uninducedData(irow, :);
        meanControl = mean(control);
        stdControl = std(control);
        ccase = inducedData(irow, :);
        prob = normcdf(mean(ccase),meanControl,stdControl);
        if prob >= 0.5
            p(irow) = 2*(1-prob);
        else p(irow) = 2*prob;
        end
        if meanControl == 0
            p(irow) = NaN;
        end
        if meanControl == 0 & mean(ccase) > smallVal
            p(irow) = 0;
        end
    end
    q = mafdr(p,'BHFDR', 'true');
    headerFinal = headerUninduced;
    dataFinal = [id, uninducedData];
    for j = 1:inducedNum
       headerFinal = [headerFinal, headerInduced{j}];
       dataFinal = [dataFinal, inducedData(:, j)];
        
    end
    headerFinal = [headerFinal, 'pVal', 'qVal'];
    dataFinal = [dataFinal, p, q];
    resultFileName = strcat(pwd, resultFile,'combinedResults/',allNames{i_file} , 'combinedResult.txt');
    fid = fopen(resultFileName, 'w');
    for j = 1:length(headerFinal)
        fprintf(fid, '%s\t', char(headerFinal{j}));
    end
    fprintf(fid, '\n');
    for i = 1:num_row
        fprintf(fid, '%d\t', dataFinal(i, 1));
        for iColumn = 2:length(headerFinal)
            fprintf(fid, '%8.10f\t', dataFinal(i, iColumn));
        end
        fprintf(fid, '\n');
        
    end
    fclose(fid);
    
end
%% End of script.










