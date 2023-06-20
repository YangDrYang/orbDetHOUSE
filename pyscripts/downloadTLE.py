# -*- coding: utf-8 -*-
"""
Script for filtering TLE based on consistency check

"""

import requests
import math
from sgp4.api import days2mdhms
from skyfield.api import load as skyfieldLoad


def tleDownload(satDict, outFile):
    """_download tle data given satellite and date information_

    Args:
        satDict (_dictionary_): _satellite information and date information_
        outFile (_string_): _*.txt file name to save tle data_
    """

    satInfo = satDict["satInfo"]
    dateInfo = satDict["dateInfo"]

    # Request URLs
    uriBase = "https://www.space-track.org"
    requestLogin = "/ajaxauth/login"
    requestCmdAction = "/basicspacedata/query"
    requestFindCustom1 = "/class/gp_history/NORAD_CAT_ID/"
    requestFindCustom2 = "/orderby/TLE_LINE1%20ASC/EPOCH/"
    requestFindCustom3 = "/format/tle"

    # API login ID
    configUsr = "yiyinfeixiong@gmail.com"
    configPwd = "Shmilyyy2008****"
    siteCred = {"identity": configUsr, "password": configPwd}

    # use requests package to drive the RESTful session with space-track.org
    with requests.Session() as session:
        # run the session in a with block to force session to close if we exit

        # need to log in first. note that we get a 200 to say the web site got the data, not that we are logged in
        resp = session.post(uriBase + requestLogin, data=siteCred)
        if resp.status_code != 200:
            print(resp)
            print("POST fail on login")

        # customSat = input("Please enter the satellite name : ")
        # catId = input("Please enter the NORAD ID of the satellite : ")
        # dateStart = input("Please enter the start date (YYYY-MM-DD format with hyphens) : ")
        # dateEnd = input("Please enter the end date (YYYY-MM-DD format with hyphens) : ")
        customSat = satInfo[0]
        catId = satInfo[1]
        dateStart = dateInfo[0]
        dateEnd = dateInfo[1]

        # Retrieves the result of the custom query
        resp = session.get(
            uriBase
            + requestCmdAction
            + requestFindCustom1
            + catId
            + requestFindCustom2
            + dateStart
            + "--"
            + dateEnd
            + requestFindCustom3
        )
        if resp.status_code != 200:
            print(resp)
            print("GET fail on request for " + customSat)

        # Assign the result line by line
        retData = resp.text.splitlines()
    # Closing the session with the API and ending the program
    session.close()

    l1 = []  # Table for line "1"
    l2 = []  # Table for line "2"
    i = 1

    # For each line of the answer, we alternate between 1 and 2
    for line in retData:
        j = i
        if i == 1:
            # If it is 1, we add the line we are reading to l1
            l1.append(line[:69])
            j = 2
        elif i == 2:
            # If it is 2, we add the line we are reading to l2
            l2.append(line[:69])
            j = 1
        i = j

    MU = (398600.4418 * 10**9) * 2 * math.pi / 86400
    k = 0
    dayList = []
    # lists for two lines, deleting extra tles and leaving one tle per day
    tleLine1List = []
    tleLine2List = []

    # print(l1)
    # print(l2)
    while k < len(l1):
        # store the days in a list to retrieve 1 data per day
        month, day, hour, minute, second = days2mdhms(
            int(l1[k][18:20]), float(l1[k][20:32])
        )

        dayList.append(day)

        if (dayList[k] != dayList[k - 1]) or k == 0:
            tleLine1List.append(l1[k])
            tleLine2List.append(l2[k])

        k = k + 1

    # reset index k to the second data
    k = 0
    # Open the outFile.txt file in write mode
    outTLEFile = open(outFile, "w")
    while k < len(tleLine1List):
        # Write the two lines to the raw file
        outTLEFile.write(tleLine1List[k] + "\n")
        outTLEFile.write(tleLine2List[k] + "\n")

        k = k + 1
    # Close the file
    outTLEFile.close()

    # Load and parse a TLE file, returning a list of Earth satellites using Skyfield API load
    dataSat = skyfieldLoad.tle_file(
        outFile,
        reload=False,
    )

    return dataSat
