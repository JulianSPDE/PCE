function w0=FUN_IC_1Ddistribution(n)

% create initial u distribution using normal distribution
mu=0; sigma=1;
p=linspace(-12,12,n);
dist=(1/(sigma*sqrt(2*pi)))*exp(-(1/2)*((p-mu)/(sigma)).^2);
dist=dist*2.5;

u0=1-(1-0.33)*dist;

% create initial v distibution using generalized normal distribution
sigma=1; var=sigma^2; alpha=sqrt(2*var);
mu=0; beta=3;

p=linspace(-14,14,n);
dist=(beta/(2*alpha*gamma(1/beta)))*exp(-(abs(p-mu)/alpha).^beta);
dist=2.5*dist;

v0=0.37*dist;

w0=[u0';v0'];

end