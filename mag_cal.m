# mag_cal: function for transforming nx3 raw measurements Mx, My, Mz
# a set of calibration constants, post-fit data, and measure of std. err.
#
# (c) 2020 David Goncalves

function [A_1, b, t, std_err] = mag_cal(s, F=1)
    # get an ellipsoid fit to data
    [M, n, d] = ellipsoid_fit(s);
    # calibration parameters A_1 and b
    M_1 = inv(M);
    b = (-M_1 * n)';
    A_1 = real(F / sqrt(n' * (M_1 * n) - d) * sqrtm(M));
    # calibrated data  and std. err.
    t = (s - b) * A_1;
    r_diff = norm(t, 2, 'rows') - 1;
    std_err = sqrt(var(r_diff)) / sqrt(rows(r_diff));
endfunction

function [M, n, d] = ellipsoid_fit (s)
    # Refer to https://teslabs.com/articles/magnetometer-calibration/
    # for derivation of magnetic calibration through an ellipsoidal fit 
    
    # D (samples)
    a1 = s(:,1) .^ 2; #a
    a2 = s(:,2) .^ 2; #b
    a3 = s(:,3) .^ 2; #c
    a4 = 2 * s(:,2) .* s(:,3); #f
    a5 = 2 * s(:,1) .* s(:,3); #g
    a6 = 2 * s(:,1) .* s(:,2); #h
    a7 = 2 * s(:,1); #p
    a8 = 2 * s(:,2); #q
    a9 = 2 * s(:,3); #r
    a10 = ones(rows(s(:,1)),1); #d
    D = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10];

    # S, S_11, S_12, S_21, S_22 (eq.11)
    S = D' * D;
    S_11 = S(1:6,1:6);
    S_12 = S(1:6,7:end);
    S_21 = S(7:end,1:6);
    S_22 = S(7:end,7:end);

    # C (Eq. 8, k=4)
    C = [ -1, 1, 1, 0, 0, 0;
          1, -1, 1, 0, 0, 0;
          1, 1, -1, 0, 0, 0
          0, 0, 0, -4, 0, 0;
          0, 0, 0, 0, -4, 0;
          0, 0, 0, 0, 0, -4;];

    # v_1 (eq. 15, solution)
    E = inv(C) * (S_11 - (S_12 * (inv(S_22) * S_21)));
    [v, l, w] = eig(E);
    E_w = diag(l);
    E_v = v;
    [m, i] = max(E_w);
    v_1 = v_1 = E_v(:,i);
    if v_1(1) < 0
      v_1 = -v_1;
    endif;

    # v_2 (eq. 13, solution)
    
    v_2 = (v_1' * (-inv(S_22) * S_21)');
      
    # quadric-form parameters
    # M = [a, h, g; h, b, f; g, f, c]
    M = [v_1(1), v_1(4), v_1(5); 
         v_1(4), v_1(2), v_1(6); 
         v_1(5), v_1(6), v_1(3)];
    # n = [p; q; r]
    n = [v_2(1); v_2(2); v_2(3)];
    d = v_2(4);

endfunction
