
function [aexp] =dbmhzcm2aexp(alpha,c0,omega0,dT)

np=-db(exp(-1));
F0=omega0/2/pi/1e6;
texp=alpha*F0*c0/1e-2/np;
aexp=exp(-dT*texp);


  if(0)

alpha=0.5;
c0=1540;
omega0=2.5e6*2*pi;
F0=omega0/2/pi/1e6;

texp=alpha*F0*c0/1e-2/np % Nepers/second
np=-db(exp(-1));

xaxis=(0:100)/100*10e-2;
taxis=xaxis/c0;
dT=taxis(2)-taxis(1);

6/(alpha*F0) %% -6 dB location in cm %%
vec1=exp(-xaxis*alpha*100*F0/np);
vec2=exp(-taxis*texp);
plot(xaxis,vec1)
grid on

hold on
plot(taxis*c0,vec2)

aexp=exp(-dT*texp)

vec3=1;
for i=2:101
  vec3(i)=vec3(i-1)*aexp;
end
plot(xaxis,vec3), hold off
  end

