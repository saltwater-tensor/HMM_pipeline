%% Sign flip correction
disp('Running sign flip correction on data...')
disp('You can pause and stop execution if you dont want to flip')
% pause(5)
options.maxlag = 10;
% options.noruns = 10;
[flips,scorepath,covmats_unflipped] = findflip(data,T,options);
data = flipdata(data,T,flips);
save('dataset_sign_flip_corrected','data','-v7.3')
save('flips_and_scorepath_and_covmats','flips','scorepath','covmats_unflipped','-v7.3')
clear options