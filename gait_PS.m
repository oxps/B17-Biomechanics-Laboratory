%% gait_PS.m -- Gait Data Analysis
%% B17 Biomechanics -- Hilary Term 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

% Load data into arrays Xmk (markers) and Xfp (forceplate) for Subject X
% You must replace 'X' below with the code for your Subject (A, B, C, D)
mk = csvread('Subject-D-markers.csv');
fp = csvread('Subject-D-forceplate.csv');

% Find the number of rows and columns in the input files
[rmk,cmk] = size(mk);    [rfp,cfp] = size(fp);

% initialize arrays of datapt number in marker file & assign to an array
% note that ':' means 'all the rows (or columns)'
datapt = zeros(rmk,1);
datapt(:,1) = mk(:,1);

% Initialize arrays for marker data (right leg only)
sac = zeros(rmk,3);      % sacrum (SACR)
asi = zeros(rmk,3);      % anterior superior iliac spine (ASIS)
thi = zeros(rmk,3);      % thigh wand marker (THI) 
kne = zeros(rmk,3);      % lateral femoral condyle (KNE)
tib = zeros(rmk,3);      % tibia wand marker (TIB)
ank = zeros(rmk,3);      % later malleolus (LMA)
hee = zeros(rmk,3);      % heel (HEE)
toe = zeros(rmk,3);      % 2nd metatarsal head (TOE)

% Assign xyz coordinates of markers to right side sacrum, asis, thigh, knee,
% ankle, heel, and toe arrays
sac(:,1:3) = mk(:,2:4);
asi(:,1:3) = mk(:,8:10);             
thi(:,1:3) = mk(:,29:31);      
kne(:,1:3) = mk(:,32:34); 
tib(:,1:3) = mk(:,35:37);
ank(:,1:3) = mk(:,38:40);      
hee(:,1:3) = mk(:,41:43);      
toe(:,1:3) = mk(:,44:46);    

%%%%%%%%%%%%% YOU NEED TO CONTINUE THE CODE FROM HERE 

% Plot yz trajectories
     figure(1)
     plot(sac(:,2),sac(:,3))
     hold on
     text(sac(rmk,2),sac(rmk,3),'SACRUM')
     plot(asi(:,2),asi(:,3))
     text(asi(rmk,2),asi(rmk,3),'ILIAC')
        
%% Question 1

    %plot the y and z trajectories of the other markers with labels
    plot(thi(:,2),thi(:,3))
    text(thi(rmk,2),thi(rmk,3),'THIGH')
    plot(kne(:,2),kne(:,3))
    text(kne(rmk,2),kne(rmk,3),'KNEE')
    plot(tib(:,2),tib(:,3))
    text(tib(rmk,2),tib(rmk,3),'TIBIA')
    plot(ank(:,2),ank(:,3))
    text(ank(rmk,2),ank(rmk,3),'ANKLE')
    plot(hee(:,2),hee(:,3))
    text(hee(rmk,2),hee(rmk,3),'HEEL')
    plot(toe(:,2),toe(:,3))
    text(toe(rmk,2),toe(rmk,3),'TOE')

    %format the graph's appearance
    axis equal
    xlabel('Y (Posterior - Anterior)')
    ylabel('Z (Inferior - Superior)')
    title('Subject D Marker Trajectories')
    legend('SACRUM','ILIAC','THIGH','KNEE','TIBIA','ANKLE','HEEL','TOE')
     
%% Question 2
     
    %calculate the vector of difference
    differences = diff(sac(:,2))/1000;

    %calculate the instantaneous velocities  
    instantaneousVelocity = differences/0.01;

    %initialise a new vector of data points
    dataPoint=(1:length(instantaneousVelocity));

    %produce plot of the instantaneous velocity against data point number
    figure(2)
    plot(dataPoint,instantaneousVelocity)
    xlabel('Data Point Number')
    ylabel('Instantaneous Velocity (metres/second)')
    title('Instantaneous Velocity of Forward Progression of Subject D')
    
    %calculate the mean velocity of the subject and plot on the graph
    averagevelocity=mean(instantaneousVelocity);
    yline(averagevelocity)

    
%% Question 3

    %calculate the distance between the knee and ankle in each dimension
    shankLengthX=abs(kne(:,1)-ank(:,1)) ;
    shankLengthY=abs(kne(:,2)-ank(:,2)) ;
    shankLengthZ=abs(kne(:,3)-ank(:,3)) ;
    
    %initialise a vector of data points starting at 1
    dataPoint=(1:length(kne));
    
    %calculate the length of the shank segments in 2D and 3D
    shankLength2D = (shankLengthY.^2 + shankLengthZ.^2).^(0.5) ;
    shankLength3D = (shankLengthX.^2 + shankLengthY.^2 + shankLengthZ.^2).^(0.5) ;
    
    %produce a plot showing the variation of the 2D length against time
    figure(3)
    subplot(2,1,1) ;
    plot(dataPoint,shankLength2D) ;
    yline(mean(shankLength2D)) ;
    xlabel('Data Point')
    ylabel('2D Shank Length (mm)')
    title('Length of 2D Shank Segment in YZ Plane')
    
    %produce a plot showing the variation of the 3D length against time
    subplot(2,1,2) ;
    plot(dataPoint,shankLength3D) ;
    yline(mean(shankLength3D)) ;
    xlabel('Data Point')
    ylabel('3D Shank Length (mm)')
    title('Length of 3D Shank Segment in XYZ Space')
   
%% Question 4(a)

    %create a new figure
    figure(4)
    hold on 

    %produce plots of the vertical trajectories of the ankle, heel and toe
    plot(dataPoint,ank(:,3))
    text(215,ank(rmk,3),'ANKLE')

    plot(dataPoint,hee(:,3))
    text(215,hee(rmk,3),'HEEL')
    
    plot(dataPoint,toe(:,3))
    text(215,toe(rmk,3),'TOE')
    
    %produce labels showing the HS, FF, HO and TO
    xline(50, '--r');
    text(45, 195, 'HS')
    
    xline(139, '--r'); 
    text(134, 195, 'HS')
    
    xline(70, '--b'); 
    text(65, 195, 'FF')
    
    xline(160, '--b'); 
    text(155, 195, 'FF')
    
    xline(83, '--g'); 
    text(78, 195, 'HO')
    
    xline(170, '--g'); 
    text(165, 195, 'HO')
    
    xline(8, '--m'); 
    text(3, 195, 'TO')
    
    xline(101, '--m'); 
    text(96, 195, 'TO')

    xline(189, '--m'); 
    text(184, 195, 'TO')
    
    %format the plots
    xlabel('Data Point')
    ylabel('Vertical (Z) Component of Marker Trajectory (mm)')
    title('Vertical Marker Trajectories of the Ankle, Heel and Toe')
    legend('ANKLE','HEEL','TOE')
    
%% Question 4(b)

    %initialise vectors to store the data to plot
    dataPoint=zeros(rfp,1);
    force = zeros(rfp,3);      % force (FOR)
    moments = zeros(rfp,3);      % moments (MOM)

    %extract the relevant data from the force plate data
    dataPoint(:,1) = fp(:,1);
    force(:,1:3) = fp(:,2:4);
    moments(:,1:3) = fp(:,5:7);
    
    %plot the ground reaction force against the data point number
    figure(5)
    hold on
    plot(dataPoint, force(:,3));
    xlabel('Data Point')
    ylabel('Vertical Component of Ground Reaction Force (N)')
    title('Vertical Component of Ground Reaction Force')
    
    %label the HS, FF, HO and TO on the plot
    xline(3735, '--r'); 
    text(3700, 245, 'HS')

    xline(3846, '--b'); 
    text(3811, 245, 'FF')

    xline(4141, '--g'); 
    text(4106, 245, 'HO')
    
    xline(4286, '--m'); 
    text(4251, 245, 'TO')
    
%% Question 6

    %calculate the vectors in the thigh and shank segments
    shankDirection=[tib(:,2), tib(:,3)] - [ank(:,2), ank(:,3)];
    thighDirection=[thi(:,2), thi(:,3)] - [kne(:,2), kne(:,3)];
    thetaValues=zeros(1,rmk);  
    
    %a loop to obtain the joint angle at each data point
    for i=1:size(shankDirection,1)
       
        %obtain the y and z position at a particular data point
        shankVector=[shankDirection(i,1); shankDirection(i,2)];
        thighVector=[thighDirection(i,1); thighDirection(i,2)];
        
        %compute the dot product of the vectors and store the value
        cosTheta = max(min(dot(shankVector,thighVector)/(norm(shankVector)*norm(thighVector)),1),-1);
        thetaValues(1,i) = acosd(cosTheta);
        
    end
    
    %plot the graph of the knee joint angle against data point number
    figure(6)
    plot(datapt, thetaValues)
    xlabel('Data Point')
    ylabel('Knee Flexion/Extension (degrees)')
    title('2D Knee Joint Angle')