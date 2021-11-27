# Geometric-Asian-Style-Call-Options-Pricing-Model
This function calculates the price of geometric Asian call options based on the closed-form solution of Kim, B. and Wee, I.S. (2014).


Function's inputs:
S0: scalar, initial price of the underlying stock;
v0: scalar, initial volatility of the stock;
theta: scalar, long run average of volatility;
sigma: scalar, the volatility of volatility;
kappa: scalar, rate of mean reversion;
rho: scalar, correlation coefficient of two brownian motions;
r: scalar, risk-free interest rate;
n: scalar, number of terms in series expansions of H and H_tilde;
T: scalar, time to maturity;
K: scalar, strike price.
Based on several testings, the price converges to be consistent as n=30. The excution speed is fast.
