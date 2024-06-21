
clear all
clc

%Load data
T1 = readtable('P9LT90_noabs_light_div.csv');
T2 = readtable('P9RT90_noabs_light_div.csv');

%% Train with P9LT, test with P9RT
Comb_mat = table2array(T1);

%Recover columns titles/names/headers
Titles = T1.Properties.VariableNames;

%Isolate step dependent parameters (4 per fly)
Comb_mat_2 = Comb_mat(:,1:5);
for i = 1:6 %6 legs
    Comb_mat_2 = [Comb_mat_2 Comb_mat(:,(32*(i-1)+10):(32*(i-1)+37))];
end
Comb_mat = Comb_mat_2;

%Pick # trials used for testing
test_trial_num = 10;

all_coeff = cell(1,1);
Tmat = Comb_mat;

%For this, remove rows of trials used of testing
newT = Tmat;
%xx = find(newT(:,2)==test_trial_num);
%newT(xx,:) = [];
Train_mat = newT;

%Fly joint features
joint_pos_1 = Train_mat; %Only 1st trial
l1 = size(joint_pos_1);
joint_pos_all = joint_pos_1; %All 1st and 2nd trial

s_joint_pos = size(joint_pos_all);
Featu_joints = [];
for i = 1:s_joint_pos(2)-5
    column_mat = joint_pos_all(:,i+5);
    pos_nan = find(isnan(column_mat));
    column_mat(pos_nan) = 0;
    Featu_joints = [Featu_joints column_mat];
end

Featu_joints_1 = Featu_joints(1:l1(1),:);

%Convert to z-score
Featu_joints_1 = zscore(Featu_joints_1);

%Remove light off (First 400 elements of each trial)
s1 = size(Featu_joints_1);
num_trial1 = s1(1)/910;

rem_ligh_off_pos_1 = [];
for i = 1:num_trial1
    xx = (910*(i-1)+1):1:(910*(i-1)+400);
    rem_ligh_off_pos_1 = [rem_ligh_off_pos_1 xx];
end
Featu_joints_1(rem_ligh_off_pos_1,:) = [];


%Ball velocities (no z-score)
ball_vel_1 = Train_mat(:,3:5);

%Remove light off
ball_vel_1(rem_ligh_off_pos_1,:) = [];

s_joint_pos = size(ball_vel_1);
all_rms = [];
t = 0.005*(1:1:510); %Time in s
i = 3; %Only z velocity

%Fit linear regression
mdl = fitlm(Featu_joints_1,ball_vel_1(:,i));
a = mdl.Coefficients;
coeff_values = a{:,1};
all_coeff{i} = coeff_values';
all_coeff_P9LT = all_coeff;





%% Train with P9RT, test with P9LT
Comb_mat = table2array(T2);

%Recover columns titles/names/headers
Titles = T1.Properties.VariableNames;

%Isolate step dependent parameters (4 per fly)
Comb_mat_2 = Comb_mat(:,1:5);
for i = 1:6 %6 legs
    Comb_mat_2 = [Comb_mat_2 Comb_mat(:,(32*(i-1)+10):(32*(i-1)+37))];
end
Comb_mat = Comb_mat_2;

%Pick # trials used for testing
test_trial_num = 10;

all_coeff = cell(1,1);
Tmat = Comb_mat;

%For this, remove rows of trials used of testing
newT = Tmat;
%xx = find(newT(:,2)==test_trial_num);
%newT(xx,:) = [];
Train_mat = newT;

%Fly joint features
joint_pos_1 = Train_mat; %Only 1st trial
l1 = size(joint_pos_1);
joint_pos_all = joint_pos_1; %All 1st and 2nd trial

s_joint_pos = size(joint_pos_all);
Featu_joints = [];
for i = 1:s_joint_pos(2)-5
    column_mat = joint_pos_all(:,i+5);
    pos_nan = find(isnan(column_mat));
    column_mat(pos_nan) = 0;
    Featu_joints = [Featu_joints column_mat];
end

Featu_joints_1 = Featu_joints(1:l1(1),:);

%Convert to z-score
Featu_joints_1 = zscore(Featu_joints_1);

%Remove light off (First 400 elements of each trial)
s1 = size(Featu_joints_1);
num_trial1 = s1(1)/910;

rem_ligh_off_pos_1 = [];
for i = 1:num_trial1
    xx = (910*(i-1)+1):1:(910*(i-1)+400);
    rem_ligh_off_pos_1 = [rem_ligh_off_pos_1 xx];
end
Featu_joints_1(rem_ligh_off_pos_1,:) = [];


%Ball velocities (no z-score)
ball_vel_1 = Train_mat(:,3:5);

%Remove light off
ball_vel_1(rem_ligh_off_pos_1,:) = [];

s_joint_pos = size(ball_vel_1);
all_rms = [];
t = 0.005*(1:1:510); %Time in s
i = 3; %Only z velocity

%Fit linear regression
mdl = fitlm(Featu_joints_1,ball_vel_1(:,i));
a = mdl.Coefficients;
coeff_values = a{:,1};
all_coeff{i} = coeff_values';
all_coeff_P9RT = all_coeff;



%% Plot coefficients
figure()

corr_coef_P9LT = all_coeff_P9LT{i};
constant_val = corr_coef_P9LT(1);
corr_coef_P9LT = corr_coef_P9LT(2:end);
subplot(1,3,1)
plot(corr_coef_P9LT,'k')
xlabel('Featured joint variables')
ylabel(['Coeff values (z velocity'])
title('Train with P9LT')
ylim([-0.6 0.6])
box off
xlim([0 length(corr_coef_P9LT)])

corr_coef_P9RT = all_coeff_P9RT{i};
constant_val = corr_coef_P9RT(1);
corr_coef_P9RT = corr_coef_P9RT(2:end);
subplot(1,3,2)
plot(corr_coef_P9RT,'k')
xlabel('Featured joint variables')
ylabel(['Coeff values (z velocity'])
title('Train with P9RT')
ylim([-0.6 0.6])
box off
xlim([0 length(corr_coef_P9RT)])


%After normalization
corr_coef_P9RT = corr_coef_P9RT/max(corr_coef_P9RT);
corr_coef_P9LT = corr_coef_P9LT/max(corr_coef_P9LT);

subplot(1,3,3)
plot(corr_coef_P9RT,'b')
hold on
plot(corr_coef_P9LT,'r')
xlabel('Featured joint variables')
ylabel(['Coeff values (z velocity'])
legend({'Train P9RT' 'Train P9LT'})
ylim([-1.2 1.2])
box off
xlim([0 length(corr_coef_P9RT)])
title('After normalization')


figure
plot(corr_coef_P9RT,'b')
hold on
plot(corr_coef_P9LT,'r')
xlabel('Featured joint variables')
ylabel(['Coeff values (z velocity'])
legend({'Train P9RT' 'Train P9LT'})
ylim([-1.5 1.2])
box off
xlim([0 length(corr_coef_P9RT)])
title('After normalization')
