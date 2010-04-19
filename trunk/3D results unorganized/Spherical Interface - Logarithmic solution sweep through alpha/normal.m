function [normal_ls] = normal(coeff)

norm_grad = sqrt (coeff(2)*coeff(2) + coeff(3)*coeff(3) + coeff(4)*coeff(4));
normal_ls(1) = -coeff(2)/norm_grad;
normal_ls(2) = -coeff(3)/norm_grad;
normal_ls(3) = -coeff(4)/norm_grad;

normal_ls = [normal_ls(1), normal_ls(2), normal_ls(3)];