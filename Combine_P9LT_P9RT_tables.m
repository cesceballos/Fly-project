
clear all
clc

%Load data
Fly_number_L = [1 2 3 4 5 7 9 10 11 12];
Fly_number_R = [1 2 3 4 5 7 8 9 10 12];

%T1 = readtable('P9LT90_noabs_light_new.csv');
%T2 = readtable('P9RT90_noabs_light_new.csv');

T1 = readtable('P9LT90_noabs_light_div.csv');
T2 = readtable('P9RT90_noabs_light_div.csv');

%Recover columns titles/names/headers
Titles = T1.Properties.VariableNames;

Tmat1 = table2array(T1);
Fly_num1 = Tmat1(:,1);
Tmat2 = table2array(T2);
Fly_num2 = Tmat2(:,1);

%Make new table combining trials per fly
Comb_mat = [];
Testing_10th_trial = [];
for i = 1:length(Fly_number_L)

    %Isolate fly data
    Fly_data1 = Tmat1(Fly_num1==Fly_number_L(i),:);
    Fly_data2 = Tmat2(Fly_num2==Fly_number_R(i),:);
    Trial_num1 = Fly_data1(:,2);
    Trial_num2 = Fly_data2(:,2);

    %Isolate trial data
    for j = 1:10 %10 trials

        Trial_data1 = Fly_data1(Trial_num1==j,:);
        Trial_data2 = Fly_data2(Trial_num2==j,:);
        Comb_mat = [Comb_mat; Trial_data1; Trial_data2];

    end
end

%save Comb_P9RT_P9LT.mat Comb_mat Titles
save Comb_P9RT_P9LT_div.mat Comb_mat Titles

