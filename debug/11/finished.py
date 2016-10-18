#!/usr/bin/python

import smtplib

sender = 'ashtonsb@umich.edu'
receivers = ['ashtonsb@umich.edu']

with open('./output/message.txt', 'r') as myfile:
    data=myfile.read()

message = """From: Ulmus <ashtonsb@ulmus.eeb.lsa.umich.edu>
To: To Person <ashtonsb@umich.edu>
Subject: CALCULATION FINISHED

Your calculation on ULMUS (ulmus.eeb.lsa.umich.edu) has finished.

""" + data

try:
   smtpObj = smtplib.SMTP('localhost')
   smtpObj.sendmail(sender, receivers, message)
   print "Successfully sent email"
except SMTPException:
   print "Error: unable to send email"
