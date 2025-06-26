function appendMatrixToCSV(filename, Hd)
    % This function appends a complex matrix Hd (MxK) to a CSV file.
    % The matrix is stored as a single row with real and imaginary parts.

    % Get the size of the matrix Hd (MxK)
    [M, K] = size(Hd);
    
    % Flatten the matrix into a 1x(M*K) row vector
    Hd_flat = reshape(Hd, 1, []); 
    
    % Separate real and imaginary parts
    Hd_flat_real = real(Hd_flat);
    Hd_flat_imag = imag(Hd_flat);

    % Combine them into alternating real and imaginary parts
    Hd_combined = [Hd_flat_real; Hd_flat_imag];
    Hd_combined = reshape(Hd_combined, 1, []); % Convert to single row
    [p, q] = size(Hd_combined);
    
    % Append the data to the CSV file
    dlmwrite(filename, Hd_combined, '-append');
    
    % Display confirmation
    fprintf('%d, %d, %d, %s\n', p*q, M, K, filename);
end
