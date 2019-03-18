%eq time constant stuff

hold on %plot all on same graph
for i=1:11 %have 11 separate things
    
    
    %create variable name
    v43 = genvarname(['t' int2str(i)]);
    eval(['xpts=' v43;]);
    
    xpts(isnan(xpts)==1)=[];
    
    %create variable name
    v43 = genvarname(['V' int2str(i)]);
    eval(['ypts=' v43;]);
    
    ypts(isnan(ypts)==1)=[];
    
    yfin(i)=ypts(end);
    
    ypts=ypts-ypts(1);
    
   [fitresult,gof]=createFit(xpts,ypts);
   
   b1(i)=fitresult.b;
   a1(i)=fitresult.a;
%   d1(i)=fitresult.d;
   c1(i)=fitresult.c;
   
end

%no fit, simple way
hold on
for i=1:3 %have 11 separate things
    
    
    %create variable name
    v43 = genvarname(['t' int2str(i)]);
    eval(['xpts=' v43;]);
       xpts(isnan(xpts)==1)=[];
    
    %create variable name
    v43 = genvarname(['V' int2str(i)]);
    eval(['ypts=' v43;]);
    ypts(isnan(ypts)==1)=[];
    
%    figure
    plot(xpts,ypts-ypts(1))

    pause
   
end