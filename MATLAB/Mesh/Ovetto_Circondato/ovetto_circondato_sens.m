function [x,y]=ovetto_circondato_sens(bs,s)

%CIRSG  Gives geometry data for the cirsg PDE model
%
%   NE=CIRSG gives the number of boundary segment
%
%   D=CIRSG(BS) gives a matrix with one column for each boundary segment
%   specified in BS.
%   Row 1 contains the start parameter value.
%   Row 2 contains the end parameter value.
%   Row 3 contains the number of the left hand region.
%   Row 4 contains the number of the right hand region.
%
%   [X,Y]=CIRSG(BS,S) gives coordinates of boundary points. BS specifies the
%   boundary segments and S the corresponding parameter values. BS may be
%   a scalar.

% Copyright 1994-2003 The MathWorks, Inc.

nbs=26;


if nargin==0,
  x=nbs; % number of boundary segments
  return
end

d=[
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 % start parameter value
  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 % end parameter value
  1 1 1 1 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 % left hand region
  2 2 2 2 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 % right hand region
];

bs1=bs(:)';

if find(bs1<1 | bs1>nbs),
  error(message('pde:cirsg:InvalidBs'))
end

if nargin==1,
  x=d(:,bs1);
  return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 && n==1,
  bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) || n~=size(s,2),
  error(message('pde:cirsg:SizeBs'));
end

if ~isempty(s),

% % boundary segment 1
% ii=find(bs==1);
% if length(ii)
% x(ii)=interp1([d(1,1),d(2,1)],[-0.70710678118654757 0],s(ii));
% y(ii)=interp1([d(1,1),d(2,1)],[0.70710678118654757 0],s(ii));
% end
% 
% % boundary segment 2
% ii=find(bs==2);
% if length(ii)
% x(ii)=interp1([d(1,2),d(2,2)],[0 -0.70710678118654757],s(ii));
% y(ii)=interp1([d(1,2),d(2,2)],[0 -0.70710678118654757],s(ii));
% end

% boundary segment 1
ii=find(bs==1);
if length(ii)
x(ii)=3.9*cos(pi/2*s(ii)+(0))+(0);
y(ii)=6.8*sin(pi/2*s(ii)+(0))+(0);
end

% boundary segment 2
ii=find(bs==2);
if length(ii)
x(ii)=3.9*cos(pi/2*s(ii)+(pi/2))+(0);
y(ii)=6.8*sin(pi/2*s(ii)+(pi/2))+(0);
end

% boundary segment 3
ii=find(bs==3);
if length(ii)
x(ii)=3.9*cos(pi/2*s(ii)+(pi))+(0);
y(ii)=3.9*sin(pi/2*s(ii)+(pi))+(0);
end

% boundary segment 4
ii=find(bs==4);
if length(ii)
x(ii)=3.9*cos(pi/2*s(ii)+(3/2*pi))+(0);
y(ii)=3.9*sin(pi/2*s(ii)+(3/2*pi))+(0);
end

% boundary segment 5
ii=find(bs==5);
if length(ii)
x(ii)=2.4*cos(pi/2*s(ii)+(0))+(0);
y(ii)=2.4*sin(pi/2*s(ii)+(0))+(0);
end

% boundary segment 6
ii=find(bs==6);
if length(ii)
x(ii)=2.4*cos(pi/2*s(ii)+(pi/2))+(0);
y(ii)=2.4*sin(pi/2*s(ii)+(pi/2))+(0);
end

% boundary segment 7
ii=find(bs==7);
if length(ii)
x(ii)=2.4*cos(pi/2*s(ii)+(pi))+(0);
y(ii)=2.4*sin(pi/2*s(ii)+(pi))+(0);
end

% boundary segment 8
ii=find(bs==8);
if length(ii)
x(ii)=2.4*cos(pi/2*s(ii)+(3/2*pi))+(0);
y(ii)=2.4*sin(pi/2*s(ii)+(3/2*pi))+(0);
end

% boundary segment 9
ii=find(bs==9);
if length(ii)
x(ii)=3.5*cos(pi*(1/2-0.1213)*(s(ii))+(pi*0.1213))+(0);
y(ii)=5.8*sin(pi*(1/2-0.1213)*(s(ii))+(pi*0.1213))+(0);
end

% boundary segment 10
ii=find(bs==10);
if length(ii)
x(ii)=3.5*cos(pi*(1/2-0.1213)*(s(ii))+(pi/2))+(0);
y(ii)=5.8*sin(pi*(1/2-0.1213)*(s(ii))+(pi/2))+(0);
end


%-3.2489,2.1571 
% boundary segment 11
ii=find(bs==11);
if length(ii)
x(ii)=3.9*cos(pi*(-1/2+0.186587647116488)*s(ii)+(pi*(1-0.186587647116488)))+(0);
y(ii)=3.9*sin(pi*(-1/2+0.186559054150777)*s(ii)+(pi*(1-0.186559054150777)))+(0);
end

% boundary segment 12
ii=find(bs==12);
if length(ii)
x(ii)=3.9*cos(pi*(-1/2+0.186587647116488)*s(ii)+(pi/2))+(0);
y(ii)=3.9*sin(pi*(-1/2+0.186559054150777)*s(ii)+(pi/2))+(0);
end

% angolo di rotazione
theta=1.5*pi/4;
x_0 = 3.5;
y_0 = -3.5;
% rho_x=15;
% rho_y=12;

rho_x = 12 ; 
rho_y = 10 ;

% boundary segment 13
ii=find(bs==13);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos(pi/2/12*s(ii)+(-pi/2/24))+(0);
y_prime=y_0+rho_y*sin(pi/2/12*s(ii)+(-pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

% boundary segment 14
ii=find(bs==14);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos((pi/3-pi/2/12)*s(ii)+(pi/2/24))+(0);
y_prime=y_0+rho_y*sin((pi/3-pi/2/12)*s(ii)+(pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

% boundary segment 15
ii=find(bs==15);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos((pi/2/12)*s(ii)+(pi/3-pi/2/24))+(0);
y_prime=y_0+rho_y*sin((pi/2/12)*s(ii)+(pi/3-pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

% boundary segment 16
ii=find(bs==16);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos((pi/2/3-pi/2/12)*s(ii)+(pi/3+pi/2/24))+(0);
y_prime=y_0+rho_y*sin((pi/2/3-pi/2/12)*s(ii)+(pi/3+pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

% boundary segment 17
ii=find(bs==17);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos((pi/2/12)*s(ii)+(pi/2-pi/2/24))+(0);
y_prime=y_0+rho_y*sin((pi/2/12)*s(ii)+(pi/2-pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

% boundary segment 18
ii=find(bs==18);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos((pi/2/6-pi/2/12)*s(ii)+(pi/2+pi/2/24))+(0);
y_prime=y_0+rho_y*sin((pi/2/6-pi/2/12)*s(ii)+(pi/2+pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

% boundary segment 19
ii=find(bs==19);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos((pi/2/12)*s(ii)+(pi/2+pi/2/6-pi/2/24))+(0);
y_prime=y_0+rho_y*sin((pi/2/12)*s(ii)+(pi/2+pi/2/6-pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

% boundary segment 20
ii=find(bs==20);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos((pi/2/6-pi/2/12)*s(ii)+(pi/2+pi/2/6+pi/2/24))+(0);
y_prime=y_0+rho_y*sin((pi/2/6-pi/2/12)*s(ii)+(pi/2+pi/2/6+pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

% boundary segment 21
ii=find(bs==21);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos((pi/2/12)*s(ii)+(pi/2+2*pi/2/6-pi/2/24))+(0);
y_prime=y_0+rho_y*sin((pi/2/12)*s(ii)+(pi/2+2*pi/2/6-pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

% boundary segment 22
ii=find(bs==22);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos((pi/2/6-pi/2/12)*s(ii)+(pi/2+2*pi/2/6+pi/2/24))+(0);
y_prime=y_0+rho_y*sin((pi/2/6-pi/2/12)*s(ii)+(pi/2+2*pi/2/6+pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

% boundary segment 23
ii=find(bs==23);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos((pi/2/12)*s(ii)+(pi/2+3*pi/2/6-pi/2/24))+(0);
y_prime=y_0+rho_y*sin((pi/2/12)*s(ii)+(pi/2+3*pi/2/6-pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

% boundary segment 24
ii=find(bs==24);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos((pi/2/2-pi/2/12)*s(ii)+(pi/2+3*pi/2/6+pi/2/24))+(0);
y_prime=y_0+rho_y*sin((pi/2/2-pi/2/12)*s(ii)+(pi/2+3*pi/2/6+pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

% boundary segment 25
ii=find(bs==25);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos((pi/2/12)*s(ii)+(pi-pi/2/24))+(0);
y_prime=y_0+rho_y*sin((pi/2/12)*s(ii)+(pi-pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

% boundary segment 26
ii=find(bs==26);
if length(ii)
s(ii)
x_prime=x_0+rho_x*cos((pi-pi/2/12)*s(ii)+(pi+pi/2/24))+(0);
y_prime=y_0+rho_y*sin((pi-pi/2/12)*s(ii)+(pi+pi/2/24))+(0);
x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
end

%%%%%%%%%%

% % boundary segment 13
% ii=find(bs==13);
% if length(ii)
%     s(ii)
%     x_prime=x_0+rho_x*cos(pi/2/6*s(ii)+(0))+(0);
%     y_prime=y_0+rho_y*sin(pi/2/6*s(ii)+(0))+(0);
% x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
% y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
% end
% 
% % boundary segment 14
% ii=find(bs==14);
% if length(ii)
%     s(ii)
%     x_prime=x_0+rho_x*cos(pi/2*5/6*s(ii)+(pi/2/6))+(0);
%     y_prime=y_0+rho_y*sin(pi/2*5/6*s(ii)+(pi/2/6))+(0);
% x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
% y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
% end
% 
% % boundary segment 15
% ii=find(bs==15);
% if length(ii)
%     x_prime=x_0+rho_x*cos(pi/2*s(ii)+(pi/2))+(0);
%     y_prime=y_0+rho_y*sin(pi/2*s(ii)+(pi/2))+(0);
% x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
% y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
% end
% 
% % boundary segment 16
% ii=find(bs==16);
% if length(ii)
%     x_prime=x_0+rho_x*cos(pi/2/6*s(ii)+(pi))+(0);
%     y_prime=y_0+rho_y*sin(pi/2/6*s(ii)+(pi))+(0);
% x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
% y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
% end
% 
% % % boundary segment 17
% ii=find(bs==17);
% if length(ii)
%     x_prime=x_0+rho_x*cos(pi/2*5/6*s(ii)+(pi+pi/2/6))+(0);
%     y_prime=y_0+rho_y*sin(pi/2*5/6*s(ii)+(pi+pi/2/6))+(0);
% x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
% y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
% end
% 
% % boundary segment 18
% ii=find(bs==18);
% if length(ii)
%     x_prime=x_0+rho_x*cos(pi/2*s(ii)+(3/2*pi))+(0);
%     y_prime=y_0+rho_y*sin(pi/2*s(ii)+(3/2*pi))+(0);
% x(ii)=x_prime*cos(theta) - y_prime*sin(theta);
% y(ii)=x_prime*sin(theta) + y_prime*cos(theta);
% end

end