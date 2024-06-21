
clear all
clc

%Load data
%load Comb_P9RT_P9LT.mat
load Comb_P9RT_P9LT_div.mat
num_flies = 10;
Fly_number_L = [1 2 3 4 5 7 9 10 11 12];
Fly_number_R = [1 2 3 4 5 7 8 9 10 12];

%Isolate step dependent parameters (4 per fly)
Comb_mat_2 = Comb_mat(:,1:5);
for i = 1:6 %6 legs
    Comb_mat_2 = [Comb_mat_2 Comb_mat(:,(32*(i-1)+6):(32*(i-1)+9))];
end
Comb_mat = Comb_mat_2;

%Pick # trials used for testing
test_trial_num = 10;

all_coeff = cell(1,1);
Tmat = Comb_mat;
Fly_num = Tmat(:,1);
Trial_num = Tmat(:,2);

%For this, remove rows of trials used of testing
newT = Tmat;
xx = find(newT(:,2)==test_trial_num);
newT(xx,:) = [];
Train_mat = newT;

%Select rows with trial used for testing
newT = Tmat;
xx = find(newT(:,2)==test_trial_num);
newT = newT(xx,:);
Test_mat = newT;

%Fly joint features
joint_pos_1 = Train_mat; %Only 1st trial
l1 = size(joint_pos_1);
joint_pos_2 = Test_mat; %Only 2nd trial
l2 = size(joint_pos_2);
joint_pos_all = [joint_pos_1; joint_pos_2]; %All 1st and 2nd trial

s_joint_pos = size(joint_pos_all);
Featu_joints = [];
featu_name = {};
featu_name_2 = {};
for i = 1:s_joint_pos(2)-5
    column_mat = joint_pos_all(:,i+5);
    pos_nan = find(isnan(column_mat));
    column_mat(pos_nan) = 0;
    Featu_joints = [Featu_joints column_mat];
end

Featu_joints_1 = Featu_joints(1:l1(1),:);
Featu_joints_2 = Featu_joints((l1(1)+1):end,:);

%Convert to z-score
Featu_joints_1 = zscore(Featu_joints_1);
Featu_joints_2 = zscore(Featu_joints_2);

%Remove light off (First 400 elements of each trial)
s1 = size(Featu_joints_1);
s2 = size(Featu_joints_2);
num_trial1 = s1(1)/910;
num_trial2 = s2(1)/910;

rem_ligh_off_pos_1 = [];
for i = 1:num_trial1
    xx = (910*(i-1)+1):1:(910*(i-1)+400);
    rem_ligh_off_pos_1 = [rem_ligh_off_pos_1 xx];
end
Featu_joints_1(rem_ligh_off_pos_1,:) = [];

rem_ligh_off_pos_2 = [];
for i = 1:num_trial2
    xx = (910*(i-1)+1):1:(910*(i-1)+400);
    rem_ligh_off_pos_2 = [rem_ligh_off_pos_2 xx];
end
Featu_joints_2(rem_ligh_off_pos_2,:) = [];


%Ball velocities (no z-score)
ball_vel_1 = Train_mat(:,3:5);
ball_vel_2 = Test_mat(:,3:5);

%Remove light off
ball_vel_1(rem_ligh_off_pos_1,:) = [];
ball_vel_2(rem_ligh_off_pos_2,:) = [];

s_joint_pos = size(ball_vel_1);
all_rms = [];
t = 0.005*(1:1:510); %Time in s
for i = 1:s_joint_pos(2)

    %Fit linear regression
    mdl = fitlm(Featu_joints_1,ball_vel_1(:,i));
    a = mdl.Coefficients;
    coeff_values = a{:,1};
    all_coeff{i} = coeff_values';
    predicted_output = Featu_joints_2*coeff_values(2:end) + coeff_values(1);

    %Calculate root mean square
    rms = sqrt(sum((ball_vel_2(:,i)-predicted_output).^2)/length(predicted_output));
    all_rms = [all_rms rms];

    figure()
    subt_pos = 0;
    for k = 1:num_flies

        subplot(2,5,k)
        plot(t,ball_vel_2((510*(k-1)+1):(510*k),i),'k')
        hold on
        plot(t,predicted_output((510*(k-1)+1):(510*k)),'r')
        ylabel('Z-score')
        ylabel('Velocity')
        xlabel('s')
        xlim([0 2.55])
        box off

        if k ==1
            %Title
            x_title = Titles{i+2};
            x_title = replace(x_title,'_',' ');
            title({x_title 'Fly #1, 10th trials'})
        end
    end
    sgtitle('P9LT')

    figure()
    subt_pos = 0;
    init_2 = 5100;
    for k = 1:num_flies

        subplot(2,5,k)
        plot(t,ball_vel_2((init_2+510*(k-1)+1):(init_2+510*k),i),'k')
        hold on
        plot(t,predicted_output((init_2+510*(k-1)+1):(init_2+510*k)),'b')
        ylabel('Z-score')
        ylabel('Velocity')
        xlabel('s')
        xlim([0 2.55])
        box off

        if k ==1
            %Title
            x_title = Titles{i+2};
            x_title = replace(x_title,'_',' ');
            title({x_title 'Fly #1, 10th trials'})
        end
    end
    sgtitle('P9RT')
end

all_rms

%Plot coefficients
figure()
veloc_vec = {'x' 'y' 'z'};
for i = 1:3

    corr_coef = all_coeff{i};
    constant_val = corr_coef(1);
    corr_coef = corr_coef(2:end);
    subplot(1,3,i)
    plot(corr_coef,'k')
    xlabel('Featured joint variables')
    ylabel(['Coeff values (' veloc_vec{i} ' velocity'])
    %ylim([-0.8 0.6])
    box off
    xlim([0 length(all_coeff{i})])
    title(['Intercept = ' num2str(constant_val)])

end


%Find correlated featured variables
corr_param = corr([Featu_joints_1; Featu_joints_2]);
corr_param(corr_param==1) = NaN;
[posx,posy] = find(corr_param>=0.9); %Highly correlated features
all_pos = [posx posy];
all_pos = unique(all_pos);

