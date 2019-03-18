% this script makes the initial support which is a rectangle with 
% half of the dimensions of the grid on which we have measured the 
% Bragg reflection and on which we will reconstruct the object.

function [ini_support] = BCDI_make_initial_support(portion)

    global X Y Z d2_bragg

    Npix_x = size(X,1);
    Npix_y = size(Y,2);
    Npix_z = size(Z,3);


    figure(1);clf;

    
    corners = [0 0 0 ; 
                Npix_x 0 0;
                0 Npix_y 0;
                Npix_x Npix_y 0;
                0 0 Npix_z ; 
                Npix_x 0 Npix_z;
                0 Npix_y Npix_z;
                Npix_x Npix_y Npix_z]*d2_bragg-[Npix_x/2 Npix_y/2 Npix_z/2]*d2_bragg;

    corners=[corners(:,1) corners(:,3) corners(:,2)];

    scatter3(corners(:,1),corners(:,3),corners(:,2),'or');
    K=convhulln(corners);
    %T=delaunayn(corners,{'Qt','Qbb','Qc','Qz'});
    %p=trisurf(K, corners(:,1), corners(:,2), corners(:,3));
    %p2=trisurf(T, corners(:,1), corners(:,2), corners(:,3));

     v1 = [corners(2,:) - corners(1,:)]; v1 = v1/norm(v1);
     v2 = [corners(3,:) - corners(1,:)]; v2 = v2/norm(v2);
     v3 = cross(v1,v2); v3 = v3/norm(v3);
     T1 = v3(1)*X + v3(2)*Y + v3(3)*Z;
     T1 = (T1>-(0.5*Npix_x-1)*d2_bragg/portion & T1<(0.5*Npix_x)*d2_bragg/portion); %two parallel lines

     Rz=[cosd(90) -sind(90) 0;
        sind(90) cosd(90) 0;
        0 0 1];
    v3_rot = (Rz*v3')';

    %  v1 = [corners(6,:) - corners(1,:)]; v1 = v1/norm(v1);
    %  v2 = [corners(4,:) - corners(1,:)]; v2 = v2/norm(v2);
    %  v3 = cross(v1,v2); 
     v3_rot = v3_rot/norm(v3_rot);
     T2 =v3_rot(1)*X + v3_rot(2)*Y + v3_rot(3)*Z;
     T2 = (T2>-(0.5*Npix_y-1)*d2_bragg/portion & T2<(0.5*Npix_y)*d2_bragg/portion); %two parallel lines

     Rx=[1 0 0;
        0 cosd(90) -sind(90);
        0 sind(90) cosd(90)];
     v3_rot = (Rx*v3')';

   
     T3 =v3_rot(1)*X + v3_rot(2)*Y + v3_rot(3)*Z;
     T3 = (T3>-(0.5*Npix_z-1)*d2_bragg/portion & T3<(0.5*Npix_z)*d2_bragg/portion); %two parallel lines

     ini_support = double(T1&T2&T3);
 
end