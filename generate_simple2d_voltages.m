Potential = SIMPLE2D_M('file1.dat');
Voltages = Potential(:, 4);
writematrix(Voltages, "Voltages.csv");