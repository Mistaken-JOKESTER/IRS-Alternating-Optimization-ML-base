function x = main
    x = 0;
    nestfunl;
    

    function nestfunl
        x = 5;
        disp("***");
        disp(x);
        disp("(((((");
    end

    x = x + 1;
end

z = main;
disp(z)