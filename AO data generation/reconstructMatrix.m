function Hd_reconstructed = reconstructMatrix(M, K, Hd_combined)
    % This function reconstructs the original complex matrix Hd (MxK) from Hd_combined

    % Split Hd_combined into real and imaginary parts
    Hd_flat_real = Hd_combined(1:2:end); % Real parts (take every other element starting from 1)
    Hd_flat_imag = Hd_combined(2:2:end); % Imaginary parts (take every other element starting from 2)

    % Combine real and imaginary parts to create a complex vector
    Hd_flat = Hd_flat_real + 1i * Hd_flat_imag;

    % Reshape the complex vector back to its original size MxK
    Hd_reconstructed = reshape(Hd_flat, M, K);
end
