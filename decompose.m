function [R, t,pts3D] = decompose(E,k,kpt1, kpt2)
        %four possible options
        [U, ~ ,V] = svd(E, 'econ');
        W = [0 -1 0;
            1  0  0;
            0  0  1];
        Re1 =   U * W * V';
        Re2 =   U * W' * V';
        te1 =   U * [0, 0, 1]';
        te1 = te1/norm(te1);
        te2 =    - te1;
        if(det(Re1) < 0) 
            Re1 = -Re1;
        end
        if(det(Re2) < 0) 
            Re2 = -Re2;
        end
        %trangulatePoint3D
        R0 = [ eye(3),  zeros(3, 1)];
        P1 = k * R0;
        P2_a = [Re1, te1];    % P2 option 1
        P2_b= [Re1, te2];    % P2 option 2
        P2_c = [Re2, te1];    % P2 option 3
        P2_d = [Re2, te2];    % P2 option  4
        P2_poten = {k*P2_a, k*P2_b, k*P2_c, k*P2_d};
        Re_poten = {Re1, Re1, Re2, Re2};
        te_poten = {te1, te2, te1, te2};
        % selecting the correct solution by triangulation of a point and choosing the solution where the point is in front of both cameras
        kpt1 = kpt1(1:2, :)';
        kpt2 = kpt2(1:2, :)';
        for i = 1 :4
                pts3D = triangulate (kpt1, kpt2,P1', P2_poten{i}')';
                p3d_1 = pts3D(:,1);
                p3d_2 = Re_poten{i} * pts3D(:,1)+ te_poten{i};
                if(p3d_1(3) > 0 && p3d_2(3) > 0)
                   R = Re_poten{i};
                   t = te_poten{i};
                   break;
                end
        end
    

end