
path = './validate_data.csv';
genData = readmatrix(path);
givenData = eph;
length(genData);
length(givenData);

posE    = genData(:, 2:4) - givenData(:, 2:4);
velE    = genData(:, 5:7) - givenData(:, 5:7);


posSE   = sum((genData(:, 2:4) - givenData(:, 2:4)).^2, 2);
velSE   = sum((genData(:, 5:7) - givenData(:, 5:7)).^2);

posRMSE = sqrt(mean(posSE))
velRMSE = sqrt(mean(velSE))

% eci2ecef = [0.158327     0.987386   0.00150628;
%             -0.987387     0.158327  -4.5087e-05;
%             -0.000283003  -0.00148014     0.999999];
        
% eph_eci = eph(:,2:4) * eci2ecef;
  

%% function fprintf('time ./validateModel.exe %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f\n', eph(1,1), eph(end,1), eph(1,2:7) 

