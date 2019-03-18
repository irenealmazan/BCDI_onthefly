global d2_bragg X Y Z

% [X Y Z] = meshgrid([-Npix/2:Npix/2-1]*d2_bragg, ...
%                     [-Npix/2:Npix/2-1]*d2_bragg,...
%                     [-depth/2:depth/2-1]*d2_bragg);

% [X Y Z] = meshgrid([-Npix/2:Npix/2-1]*d2_bragg, ...
%                     [-Npix/2:Npix/2-1]*d2_bragg,...
%                     [-depth/2:depth/2-1]*th_pixel_size);


figure(1);clf;
hold on;

corepix = corewidth/d2_bragg;
lenpix = round(length/d2_bragg);

edge = corewidth;
facet_spacing = 2*sqrt((edge/2)^2 - (edge/4)^2);

corners = [ edge/2 0 length/2; 
            edge/2*cosd(60) edge/2*sind(60) length/2;
            -edge/2*cosd(60) edge/2*sind(60) length/2;
            -edge/2 0 length/2;
            edge/2*cosd(60) -edge/2*sind(60) length/2;
            -edge/2*cosd(60) -edge/2*sind(60) length/2;
            edge/2 0 -length/2;
            edge/2*cosd(60) edge/2*sind(60) -length/2;
            -edge/2*cosd(60) edge/2*sind(60) -length/2;
            -edge/2 0 -length/2;
            edge/2*cosd(60) -edge/2*sind(60) -length/2;
            -edge/2*cosd(60) -edge/2*sind(60) -length/2];
corners=[corners(:,1) corners(:,3) corners(:,2)];
corners_o = corners;

Rz=[cosd(phi) -sind(phi) 0;
    sind(phi) cosd(phi) 0;
    0 0 1];
corners = (Rz*corners')';

Ry=[cosd(th) 0 sind(th);
    0 1 0;
    -sind(th) 0 cosd(th)];
corners = (Ry*corners')';

Ry=[cosd(-del) 0 sind(-del);
    0 1 0;
    -sind(-del) 0 cosd(-del)];
corners = (Ry*corners')';

Rx=[1 0 0;
    0 cosd(gam) -sind(gam);
    0 sind(gam) cosd(gam)];
corners = (Rx*corners')';

K=convhulln(corners);
T=delaunayn(corners,{'Qt','Qbb','Qc','Qz'});
p=trisurf(K, corners(:,1), corners(:,2), corners(:,3));
set(p,'FaceColor','red','EdgeColor','black');
alpha(.3);
%set(gca, 'Projection', 'perspective');
axis equal
xlabel('x (microns)'); ylabel('y'); zlabel('z');
hold on; 
quiver3(0,0,0, 0, 0, .2);
drawnow

hold off;

%%

%make parallel planes parallel to each pair of facets
v1 = [corners(2,:) - corners(8,:)]; v1 = v1/norm(v1);
v2 = [corners(1,:) - corners(8,:)]; v2 = v2/norm(v2);
v3 = cross(v1,v2); v3 = v3/norm(v3);
T1 = v3(1)*X + v3(2)*Y + v3(3)*Z;
T1 = (T1>-facet_spacing/2 & T1<facet_spacing/2); %two parallel lines

v1 = [corners(2,:) - corners(8,:)]; v1 = v1/norm(v1);
v2 = [corners(3,:) - corners(8,:)]; v2 = v2/norm(v2);
v3 = cross(v1,v2); v3 = v3/norm(v3);
T2 = v3(1)*X + v3(2)*Y + v3(3)*Z;
T2 = (T2>-facet_spacing/2 & T2<facet_spacing/2); %two parallel lines

v1 = [corners(4,:) - corners(9,:)]; v1 = v1/norm(v1);
v2 = [corners(3,:) - corners(9,:)]; v2 = v2/norm(v2);
v3 = cross(v1,v2); v3 = v3/norm(v3);
T3 = v3(1)*X + v3(2)*Y + v3(3)*Z;
T3 = (T3>-facet_spacing/2 & T3<facet_spacing/2); %two parallel lines

v3 = [0 1 0];
T4 = v3(1)*X + v3(2)*Y + v3(3)*Z;
T4 = (T4>-facet_spacing/2 & T4<facet_spacing/2); %two parallel lines

NW = double(T1&T2&T3&T4);
NWbase = NW;

if (addNWstrain)
    dispfield = 8*(sin(X) + cos(Z+Y).^2 + tan(Y));
    NW = NW.*exp(i*dispfield);
end

NWsfzb = NW;

if (addNWsf)
    
    v1 = [corners(1,:) - corners(7,:)]; v1 = v1/norm(v1);
    Tsf = v1(1)*X + v1(2)*Y + v1(3)*Z;
    sfheight = .004;
    sfboundaries = [-length/2:sfheight:length/2];
    sfboundaries = sfboundaries + (rand(size(sfboundaries))-.5) * .75*sfheight;
    %sfphases = rand(size(sfboundaries)) * 2*pi - pi;
    sfphases = round(rand(size(sfboundaries))*2-1) *2/3*pi;
    sfphases2 = sign(-.5+rand(size(sfboundaries))) *2/3*pi;
    sfamp = round(rand(size(sfboundaries))*2);
    sfamp(find(sfamp==2))=1;
    
    hold on; h1=[]; 
    for ii=2:numel(sfboundaries)
        
        ind = find(Tsf>sfboundaries(ii-1) & Tsf<=sfboundaries(ii) & abs(NWbase)==1);
        display(['found ' num2str(numel(ind)) ' in SF boundary']);
        
        NWsfzb(ind) = NW(ind)*exp(i*sfphases(ii)) * sfamp(ii);
        %NWsfzb(ind) = NW(ind)*exp(i*(sfphases2(ii)+sfphases2(ii-1))) * sfamp(ii);
        NW(ind) = NW(ind)*exp(i*sfphases(ii));
        
        delete(h1);
        temp = zeros(size(NW)); temp(ind) = sfamp(ii);
        h1 = di(temp, -.5, 'g', X,Y,Z);
        pause(.1);
        
    end
    hold off;
    
end
