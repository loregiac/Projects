clear all
close all

dt = 4;
imgap = 1600;
nim = 1080;

for i = 1 : nim
    temp = dlmread(['temperatures/temperature',int2str(imgap*(i)),'.dat']);   
    figure(1)
    plot(temp(:,1),temp(:,2),'r')
    hold on
    plot(temp(:,1),temp(:,3),'b')
    legend('Tf','Ts')
    title(int2str(imgap*(i)*dt))
    grid on
    ylim([200,900])
    title(['t = ',int2str(imgap*(i)*dt),'s'])
    hold off    
end
