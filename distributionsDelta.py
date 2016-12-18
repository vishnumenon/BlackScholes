from yahoo_finance import Share
from pprint import pprint
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import statistics as stat
from datetime import datetime


def getPriceOnDay(share, date):
	data = share.get_historical(date, date)
	return (float(data[0]['High']) + float(data[0]['Low']))/2.0

def formatDate(date):
	date = str(date)
	date = date[0:4] + "-" + date[4:6] + "-" + date[6:8]
	return date

def normalizeFreqDictionary(dict):
	total_sum = sum(dict.values())
	factor = 1.0/total_sum
	return {key:value*factor for key,value in dict.iteritems()}

def normalizeDictionary(dict, counterdict):
	for key,value in dict.iteritems():
		fac = counterdict[key]
		dict[key] = value/fac
	total_sum = sum(dict.values())
	factor = 1.0/total_sum
	return {key:value*factor for key,value in dict.iteritems()}

def reconstructDeltaValues(deltas, last_index, next_index, old_value, new_value):
	factor = float((old_value - new_value)/(next_index - last_index))
	curr_newval = old_value - factor
	for i in range(last_index + 1, next_index):
		deltas[i] = curr_newval
		curr_newval = curr_newval - factor
	return deltas

def days_between(d1, d2):
    d1 = datetime.strptime(str(d1), "%Y%m%d")
    d2 = datetime.strptime(str(d2), "%Y%m%d")
    return abs((d2 - d1).days)

first_date = "2015-09-01"
end_date = "2016-03-01"
ticker = 'AAPL'

stock_object = Share(ticker)
df = pd.read_csv(os.path.abspath("AAPL.csv"))
df['strike_price']*=0.001

expiry_dates = sorted(df['exdate'].drop_duplicates().tolist())
dates = sorted(df['date'].drop_duplicates().tolist())

dist_dict = dict()
counter_dict = dict()
realdist_dict = dict()
days_left = set()
counter = 0

for exdate in expiry_dates:
	if counter > 10:
		break
	counter += 1
	for date in dates:
		# Now we have the section of data for a particular date and an expiration
		if date >= exdate:
			break
		snip = df[(df.exdate == exdate) & (df.date == date) & (df.cp_flag == "C")]
		snip = snip.drop_duplicates(subset = 'strike_price')
		if len(snip) == 0:
			continue
		snip = snip.sort('strike_price', ascending=1)
		# Now we take the range of delta values and construct a distribution

		daysToExpiration = days_between(date, exdate)
		days_left.add(daysToExpiration)
		current_price = getPriceOnDay(stock_object, formatDate(date))

		delta_values = snip['delta'].tolist()
		strike_values = snip['strike_price'].tolist()
		original_strikes = snip['strike_price'].tolist()

		# Mean center strike values
		strike_values = [round(sp - current_price) for sp in strike_values]
		print daysToExpiration
		
		# Some parts of the delta column are incomplete in the WRDS data sheets
		# Need to make an assumption on the delta values, interpolate and fit to 1 0 spread

		last_index = -1
		last_value = 1
		for index, d in enumerate(delta_values):
			if np.isnan(float(d)):
				continue
			else:
				if index != last_index + 1:
					delta_values = reconstructDeltaValues(delta_values, last_index, index, last_value, d)
				last_index = index
				last_value = delta_values[index]
		# finish cleaning up tail values
		if len(delta_values) > last_index + 1:
			factor = float(last_value/(len(delta_values) - 1 - last_index))
			curr_val = last_value - factor
			for i in range(last_index+1, len(delta_values)):
				delta_values[i] = curr_val
				curr_val = curr_val - factor


		# Now, delta value are cleaned up and ready to be bucketed
		remaining_p = 1
		for index, d in enumerate(delta_values):
			delta_values[index] = remaining_p - d
			remaining_p = d

		# remove negative delta values (assumption)

		for i, d in enumerate(delta_values):
			if d < 0:
				delta_values[i] = 0

		future = getPriceOnDay(stock_object, formatDate(exdate))

		just_itm = round(max([contract for contract in original_strikes if contract < future]) - current_price)
		print "CONTRACT IS", just_itm

		if dist_dict.has_key(daysToExpiration):
			dist = dist_dict[daysToExpiration]
			realdist = realdist_dict[daysToExpiration]
			counterdist = counter_dict[daysToExpiration]
			for index, sp in enumerate(strike_values):
				if not dist.has_key(sp):
					dist[sp] = delta_values[index]
					counterdist[sp] = 1
				else:
					dist[sp] += delta_values[index]
					counterdist[sp] += 1
			dist_dict[daysToExpiration] = dist
			counter_dict[daysToExpiration] = counterdist
			if not realdist.has_key(just_itm):
				realdist[just_itm] = 1
			else:
				realdist[just_itm] += 1
			realdist_dict[daysToExpiration] = realdist
		else:
			dist = dict()
			realdist = dict()
			counterdist = dict()
			for index, sp in enumerate(strike_values):
				dist[sp] = delta_values[index]
				counterdist[sp] = 1
			dist_dict[daysToExpiration] = dist
			realdist[just_itm] = 1
			realdist_dict[daysToExpiration] = realdist
			counter_dict[daysToExpiration] = counterdist

# plotting probabilities

days_left =  list(days_left)
number_dists = len(days_left)
print number_dists

targets = [days_left[0], days_left[5], days_left[number_dists -1]]

for index, element in enumerate(targets):
	dist_dict[element] = normalizeDictionary(dist_dict[element], counter_dict[element])
	plt.subplot(len(targets),1,index+1)
	plt.axis([-60, 60, 0, .3])
	plt.title(element)
	plt.bar(dist_dict[element].keys(), dist_dict[element].values())
	counter = counter + 1
plt.show()

for index, element in enumerate(targets):
	pprint(realdist_dict[element])
	realdist_dict[element] = normalizeFreqDictionary(realdist_dict[element])
	plt.subplot(len(targets),1,index+1)
	plt.axis([-60, 60, 0, .3])
	plt.title(element)
	plt.bar(realdist_dict[element].keys(), realdist_dict[element].values())
	counter = counter + 1
plt.show()

# for index, element in enumerate(targets):
#realdist_dict[element] = normalizeDictionary(realdist_dict[element])
# 	plt.subplot(len(targets),1,index+1)
# 	plt.axis([-60, 60, 0, .3])
# 	plt.title(element)
	# plt.bar(realdist_dict[element].keys(), realdist_dict[element].values())
# 	counter = counter + 1
# plt.show()










		


