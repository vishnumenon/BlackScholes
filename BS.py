from scipy.stats import norm
from math import *
from datetime import datetime, date, timedelta
import requests
import json
import re

# T: time to maturity
# S: Spot price
# K: Strike price
# r: risk free rate
# sigma: volatility of returns of underlying asset
# d: yield on asset
def blackScholes(T, S, K, r, sigma, d):
    d1 = (log(float(S) / K) + (r-d)*T + sigma*sigma*(T/2.0))/float(sigma*sqrt(T))
    d2 = d1 - sigma*sqrt(T)
    C = S*exp(-1*d*T)*norm.cdf(d1) - K*exp(-1*r*T)*norm.cdf(d2)
    P = K*exp(-1*r*T)*norm.cdf(-1*d2) - S*exp(-1*d*T)*norm.cdf(-1*d1)
    return {
        "call": {
            "price": C,
            "delta": norm.cdf(d1),
            "gamma": norm.pdf(d1) / (S*sigma*sqrt(T)),
            "vega": S*norm.pdf(d1)*sqrt(T),
            "theta": (-1.0*S*norm.pdf(d1)*sigma / (2*sqrt(T))) - r*K*exp(-1*r*T)*norm.cdf(d2),
            "rho": K*T*exp(-1*r*T)*norm.cdf(d2)
        },
        "put": {
            "price": P,
            "delta": norm.cdf(d1) - 1,
            "gamma": norm.pdf(d1) / (S*sigma*sqrt(T)),
            "vega": S*norm.pdf(d1)*sqrt(T),
            "theta": (-1.0*S*norm.pdf(d1)*sigma / (2*sqrt(T))) + r*K*exp(-1*r*T)*norm.cdf(-1*d2),
            "rho": -1*K*T*exp(-1*r*T)*norm.cdf(-1*d2)
        }
    }

# TBill as Proxy
billRatesResponse = requests.get("https://www.quandl.com/api/v3/datasets/USTREASURY/BILLRATES.json?api_key=xfcug5J4ks3ZNtL5DGuo").json()
currentData = billRatesResponse["dataset"]["data"][0]
riskFreeRates = {4: currentData[2], 13: currentData[4], 26: currentData[6], 52: currentData[8]}
# Consider using Swaps as proxy? or LIBOR?

def getRiskFreeRate(years):
    maturityBucket = min(riskFreeRates.keys(), key=lambda weeks: abs(years*52 - weeks))
    return riskFreeRates[maturityBucket]

# initial estimate from Brenner-Subrahmanyam (1988)
# Consider using Minqiang Li (2006)?
# then iterating w/ newton's method
def getImpliedVolatility(T, S, C, cStrike, r, yieldValue):
    sigma = -100
    sigmaNew = sqrt(2.0*pi / T) * (float(C) / S)
    while abs(sigmaNew - sigma) > 0.000001:
        sigma = sigmaNew
        bs = blackScholes(T, S, cStrike, r, sigma, yieldValue)
        sigmaNew = sigma - ((bs["call"]["price"] - C) / bs["call"]["vega"])
    return sigmaNew

def getSpotPrice(ticker):
    response = requests.get("https://www.google.com/finance/info?q=" + ticker)
    response = json.loads(response.content[4:-1])
    return float(response[0]["l"])

def getAtMoneyCall(ticker, expirationDate, spotPrice):
    response = requests.get("http://www.google.com/finance/option_chain?q=" + ticker + "&expd=" + str(expirationDate.day) + "&expm=" + str(expirationDate.month) + "&expy=" + str(expirationDate.year) + "&output=json")
    response = response.text
    response = re.sub('(\w+):', r'"\1":', response)
    response = json.loads(response)
    pricedCalls = filter(lambda call: call["p"] != "-", response['calls'])
    return min(pricedCalls, key=lambda call: abs(float(call["strike"].replace(",", "")) - spotPrice))

def ceilExpirationDate(expirationDate):
    while not (14 < expirationDate.day < 22 and expirationDate.weekday() == 4):
        expirationDate = expirationDate + timedelta(days=1)
    return expirationDate

def priceOption(ticker, strikePrice, expiration, yieldValue = 0):
    today = date.today()
    expirationDate = datetime.strptime(expiration, '%b %d %Y')
    fixedExpirationDate = ceilExpirationDate(expirationDate).date()
    S = getSpotPrice(ticker)
    atMoneyCall = getAtMoneyCall(ticker, fixedExpirationDate, S)
    C = float(atMoneyCall["p"].replace(",", ""));
    cStrike = float(atMoneyCall["strike"].replace(",", ""));
    T = (fixedExpirationDate - today).days / 365.0
    r = getRiskFreeRate(T);
    sigma = getImpliedVolatility(T, S, C, cStrike, r, yieldValue);
    return blackScholes(T, S, strikePrice, r, sigma, yieldValue)

print(priceOption("AMZN", 730, "Oct 20 2016"));
