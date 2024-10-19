for i=2:11
    s = ["Test_matrix_" num2str(i) ".csv"];
    x = dlmread(s, " ");
    chol(x)
    disp([num2str(i) ": " num2str(det(x))])
end
