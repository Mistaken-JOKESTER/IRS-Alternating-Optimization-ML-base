iterations = 20;

for i=1:1000
    disp(['Counter is: ', num2str(i)]);
    run converge_A2_perfect.m;
end