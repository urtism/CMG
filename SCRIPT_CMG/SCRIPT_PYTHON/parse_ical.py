import argparse
import re

	
calendar = {
	'id':'',
	'version':'',
	'scale':'',
	'method':'',
	'name':'',
	'timezone':'',
	'timezone_arr':[],
	'event_arr':[]
}

event = {
	'DTSTART':'',
	'DTSEND':'',
	'DTSTAMP':'',
	'UID':'',
	'CREATED':'',
	'DESCRIPTION':'',
	'LAST-MODIFIED':'',
	'LOCATION':'',
	'SEQUENCE':'',
	'STATUS':'',
	'SUMMARY':'',
	'TRANSP':''
}

timezone = {

	'TZID':'',
	'X-LIC-LOCATION':'',
	'DAYLIGHT':{
		'TZOFFSETFROM':'',
		'TZOFFSETTO':'',
		'TZNAME':'',
		'DTSTART':'',
		'RRULE':{
			'FREQ':'',
			'BYMONTH':'',
			'BYDAY':''
		}
	},

	'STANDARD':{
		'TZOFFSETFROM':'',
		'TZOFFSETTO':'',
		'TZNAME':'',
		'DTSTART':'',
		'RRULE' : {
			'FREQ':'',
			'BYMONTH':'',
			'BYDAY':''
		}
	}
}

if __name__ == '__main__':

	parser = argparse.ArgumentParser('Parse ics file in tab delimited')
	parser.add_argument('-c', '--calendar', help="File ics da parsare")
	parser.add_argument('-o', '--out', help="file csv parsato")

	global opts
	opts = parser.parse_args()
	calendarlines=[]

	calendar = {
		'id':'',
		'version':'',
		'scale':'',
		'method':'',
		'name':'',
		'timezone':'',
		'timezone_arr':[],
		'event_arr':[]
	}

	out = open(opts.out,'w')
	
	for a in open(opts.calendar,'r').readlines():
		calendarlines += [a.strip('\r\n')] 

	#for line in calendarlines:
	#	calendar_text=calendarlines[calendarlines.index('BEGIN:VCALENDAR'):calendarlines.index('END:VCALENDAR')]
		#for ev in calendar_text.index('BEGIN:VEVENT'):
		#print calendar_text.index('BEGIN:VEVENT')

	

	start_event_indexes = [index for index in range(len(calendarlines)) if calendarlines[index] == 'BEGIN:VEVENT']
	end_event_indexes = [index for index in range(len(calendarlines)) if calendarlines[index] == 'END:VEVENT']

	#print len(start_event_indexes)
	#print len(end_event_indexes)

	event_arr=[]

	for istart in start_event_indexes:

		iend = end_event_indexes[start_event_indexes.index(istart)]

		event_arr+=[calendarlines[istart+1:iend]]

	for ev in event_arr:

		event = {
			'DTSTART':'',
			'DTEND':'',
			'DEND':'',
			'TEND':'',
			'TSTART':'',
			'DTSEND':'',
			'TSEND':'',
			'DTSTAMP':'',
			'TSTAMP':'',
			'UID':'',
			'CREATED':'',
			'DESCRIPTION':'',
			'LAST-MODIFIED':'',
			'LOCATION':'',
			'SEQUENCE':'',
			'STATUS':'',
			'SUMMARY':'',
			'TRANSP':''
		}

		for tag in ev:
			if tag.startswith('DTSTART'):
				try:
					val=tag.split(':')[1]
					event['DTSTART']=val.split('T')[0]
					try:
						event['TSTART']=val.split('T')[1]
					except:
						event['TSTART']=''

				except:
					event['DTSTART']=''
					event['TSTART']=''

			if tag.startswith('DTEND'):
				try:
					val=tag.split(':')[1]
					event['DTEND']=val.split('T')[0]
					try:
						event['TEND']=val.split('T')[1]
					except:
						event['TEND']=''
				except:
					event['DTEND']=''
					event['DEND']=''


			if tag.startswith('DTSTAMP'):

				try:
					val=tag.split(':')[1]
					event['DTSTAMP']=val.split('T')[0]
					try:
						event['TSTAMP']=val.split('T')[1]
					except:
						event['TSTAMP']=''

				except:
					event['DTSTAMP']=''
					event['TSTAMP']=''

			if tag.startswith('UID'):
				try:
					event['UID']=tag.split(':')[1]
				except:
					event['UID']=''

			if tag.startswith('CREATED'):
				try:
					event['CREATED']=tag.split(':')[1]
				except:
					event['CREATED']=''

			if tag.startswith('DESCRIPTION'):
				try:
					event['DESCRIPTION']=tag.split(':')[1]
				except:
					event['DESCRIPTION']=''

			if tag.startswith('LAST-MODIFIED'):
				try:
					event['LAST-MODIFIED']=tag.split(':')[1]
				except:
					event['LAST-MODIFIED']=''

			if tag.startswith('LOCATION'):
				try:
					event['LOCATION']=tag.split(':')[1]
				except:
					event['LOCATION']=''

			if tag.startswith('SEQUENCE'):
				try:
					event['SEQUENCE']=tag.split(':')[1]
				except:
					event['SEQUENCE']=''

			if tag.startswith('STATUS'):
				try:
					event['STATUS']=tag.split(':')[1]
				except:
					event['STATUS']=''

			if tag.startswith('SUMMARY'):
				try:
					event['SUMMARY']=tag.split(':')[1]
				except:
					event['SUMMARY']=''

			if tag.startswith('TRANSP'):
				try:
					event['TRANSP']=tag.split(':')[1]
				except:
					event['TRANSP']=''

		calendar['event_arr']+=[event]
			
	out.write('\t'.join(['DTSTART','TSTART','DTEND','TEND','DTSTAMP','TSTAMP','UID','CREATED','DESCRIPTION','LAST-MODIFIED','LOCATION','SEQUENCE','STATUS','SUMMARY'])+'\n')

	for event in calendar['event_arr']:
		try:
			out.write('\t'.join([event['DTSTART'],event['TSTART'],event['DTEND'],event['TEND'],event['DTSTAMP'],event['TSTAMP'],event['UID'],event['CREATED'],event['DESCRIPTION'],
				event['LAST-MODIFIED'],event['LOCATION'],event['SEQUENCE'],event['STATUS'],event['SUMMARY']])+'\n')
		except:
			print event
	#print event_arr