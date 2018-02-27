clear all
close all
dt = 4;
imgap = 1600;
nim = 1080;

for i = 1 : nim
    temp = dlmread(['temperatures/temperature',int2str(imgap*(i)),'.dat']);
    N = size(temp,1);
    Tf = repmat(temp(:,2),1,100);
    Ts = repmat(temp(:,3),1,100);
    [x,y] = meshgrid(1:100,linspace(temp(N,1),0,N)); 
    figure(1)
    subplot(1,2,1) 
    surf(x,y,Tf,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
    hold on
    colormap jet
    colorbar
    caxis([293,873])
    axis([0,100,0,temp(N,1),0,1000])
    view([0,0,1])
    title(['Fluid phase -- t = ',int2str(imgap*(i)*dt),'s'])
    hold off
    subplot(1,2,2) 
    surf(x,y,Ts,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
    hold on
    colormap jet
    colorbar
    caxis([293,873])
    axis([0,100,0,temp(N,1),0,1000])
    view([0,0,1])
    title(['Solid phase -- t = ',int2str(imgap*(i)*dt),'s'])
    hold off
end
