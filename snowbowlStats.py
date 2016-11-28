GlowScript 2.2 VPython
#data from my program azSnowBowlStats.py which pulls data from arizonasnowbowl's website
#and forecast.weather.gov (using BeautifulSoup) to get info related to the snow season, including snow depth and
#temperature, data in above program is managed using Pandas DataFrame
#olddata is 2015-2016 season
#data is the 2016-2017 season
#data[0] is the 0th day's entries from data[0]
#cols are the column headers for these entries
#olddata is missing 3 numbers: LowTemp,HighTemp,PrecipProb, and Temp becomes AzSbTemp
#future data will be stored in another dimension of data if its structure stays the same
#The task of plotting this was made most difficult by glowscript not slicing with 2+dimensions

#I CAN'T PLOT VS A STRING/DATE, ATTEMPTING TO WORK AROUND THIS WITH CONCATENATION
#LEADS TO JUMPS AT THE TURN OF THE YEAR. THE BEST OPTION IS TO PLOT VS 
#DECIMALS OF A MONTH, RELATIVE TO JAN 1ST (X=0)

#2016-2017 season
cols=['Date','Weekday','MinDepth','MaxDepth','NewBaseSnow','NewPeakSnow','AzSbTemp','LowTemp','HighTemp','PrecipProb','OpenLifts','Lifts','OpenTrails','Trails','SnowCond','RoadCond','WeatherCond','SnowHazard','RoadHazard']
data = [['2016-11-25','Friday',17,17,0,0,29,float('NaN'),float('NaN'),float('NaN'),1,5,7,47,'POWDER/PACKED POWDER/GROOMED','CLEAR','SUNNY AND BEAUTIFUL WEATHER FORECAST THROUGH SATURDAY! MORE SNOW ON THE WAY FOR SUNDAY AND MONDAY!',3,0],
['2016-11-26','Saturday',17,17,0,0,32,float('NaN'),float('NaN'),float('NaN'),1,5,7,47,'POWDER/PACKED POWDER/GROOMED','CLEAR','SUNNY AND BEAUTIFUL WEATHER FORECAST THROUGH SATURDAY! MORE SNOW ON THE WAY FOR SUNDAY AND MONDAY!!',3,0],
['2016-11-27','Sunday',18,18,1,1,15,17,35,100,1,5,7,47,'POWDER/ PACKED POWDER','CHAINS/ 4WD','1 INCH OF NEW SNOW OVERNIGHT AND MORE SNOW ON THE WAY THROUGH TODAY AND MONDAY.',2,4],
['2016-11-28','Monday',22,22,6,6,22,9,24,70,1,5,7,47,'POWDER/ PACKED POWDER','CHAINS/ 4WD','6 INCHES OF SNOW IN THE LAST 24 HOURS. 7 INCH STORM TOTAL.',2,4]
]

#2015-2016 season
oldcols=['Date','Weekday','MinDepth','MaxDepth','NewBaseSnow','NewPeakSnow','Temp','OpenLifts','Lifts','OpenTrails','Trails','SnowCond','RoadCond','WeatherCond','SnowHazard','RoadHazard']
olddata = [['2015-11-22','Sunday',24,33,0,0,48,float('NaN'),float('NaN'),float('NaN'),float('NaN'),'Groomed variable conditions','Clear','Sunny skies  with high of 48 with winds out of the NE 9 mph.',4,0],
['2015-11-23','Monday',24,33,0,0,48,float('NaN'),float('NaN'),float('NaN'),float('NaN'),'Groomed variable conditions','Clear','Sunny skies  with high of 48 with light winds out of the SW 8-13 mph.',4,0],
['2015-11-24','Tuesday',24,33,0,0,44,float('NaN'),float('NaN'),float('NaN'),float('NaN'),'Groomed variable conditions','Clear','Sunny skies  with high of 44 with moderate winds out of the SW 12-20 mph.',4,0],
['2015-11-25','Wednesday',24,33,0,0,24,float('NaN'),float('NaN'),float('NaN'),float('NaN'),'Groomed variable conditions','Clear','20% chance of am snow and 40% chance tonight. Winds out of the SW 20-32 mph.',4,0],
['2015-11-26','Thursday',24,33,0,0,29,float('NaN'),float('NaN'),float('NaN'),float('NaN'),'Groomed variable conditions','Clear','Sunny with a high near 29 and calm SW winds 7-14 mph.',4,0],
['2015-11-27','Friday',24,30,0,0,25,float('NaN'),float('NaN'),float('NaN'),float('NaN'),'Groomed variable conditions','Clear','Increasing clouds 40% chance of snow in the afternoon/evening.',4,0],
['2015-11-28','Saturday',24,30,1,2,26,float('NaN'),float('NaN'),float('NaN'),float('NaN'),'Groomed variable conditions','4WD or chains recommended','Sunny a slight chance of showers in the evening',4,3],
['2015-11-29','Sunday',24,30,0,0,27,3,5,24,47,'Groomed variable conditions','4WD or chains recommended','high near 27 with a 20% chance of snow calm SW winds 5-13 mph',4,3],
['2015-11-30','Monday',20,28,0,0,33,3,5,24,47,'Groomed variable conditions','Clear','high near 33 with calm-moderate winds out of the west 8-13 mph with possible am gusts.',4,0],
['2015-12-01','Tuesday',20,28,0,0,35,3,5,23,47,'Groomed variable conditions','Clear','high of 35 with calm NE winds 0-8 mph. Great Bluebird Day!',4,0],
['2015-12-02','Wednesday',20,28,0,0,42,3,5,23,47,'Groomed variable conditions','Clear','high of 42 with calm E winds 0-5 mph. Another Great Bluebird Day!',4,0],
['2015-12-03','Thursday',20,28,0,0,46,3,5,23,47,'Groomed variable conditions','Clear','high of 46 with calm SE  winds 5-14 mph. Another Great Bluebird Day!',4,0],
['2015-12-04','Friday',20,28,0,0,41,3,5,23,47,'Groomed variable conditions','Clear','high of 41 with calm-moderate SW winds 12-20 mph.',4,0],
['2015-12-05','Saturday',20,28,0,0,42,3,5,23,47,'Groomed variable conditions','Clear','high of 42 with north-west wind 5-13 mph',4,0],
['2015-12-06','Sunday',20,28,0,0,49,3,5,23,47,'Groomed variable conditions','Clear','For Sunsay high of 49 with southwest wind 9-11 mph',4,0],
['2015-12-07','Monday',20,28,0,0,50,3,5,23,47,'Groomed variable conditions','Clear','high of 50 with southwest wind 5-11 mph',4,0],
['2015-12-08','Tuesday',20,28,0,0,50,3,5,22,47,'Groomed variable conditions','Clear','high of 50 with northwest wind around 6 mph',4,0],
['2015-12-09','Wednesday',20,28,0,0,52,3,5,22,47,'Groomed variable conditions','Clear','high of 55 with north wind 5-10 mph',4,0],
['2015-12-10','Thursday',20,28,0,0,51,3,5,22,47,'Groomed variable conditions','Clear','high of 51 with west wind 10-20 mph',4,0],
['2015-12-11','Friday',20,28,0,0,34,3,5,22,47,'Groomed variable conditions','Clear','high of 34 and dropping to mid 20\'s by 5 pm with SW wind 12-28 mph breezy at times.',4,0],
['2015-12-12','Saturday',20,28,1,2,25,3,5,22,47,'Fresh powder and groomed variable conditions','4WD or chains REQUIRED','2\" new snow. High of 25 with wind 6-14 mph. Chance of precipitation 90% - another 1-3\" possible',2,4],
['2015-12-13','Sunday',30,40,10,12,38,3,5,26,47,'Fresh powder and groomed variable conditions','4WD or chains RECOMMENDED','Agassiz currently on Wind Hold.',2,3],
['2015-12-14','Monday',30,36,1,3,20,3,5,26,47,'Fresh powder and groomed variable conditions','4WD or chains REQUIRED','3\" of fresh snow another 5-9\" possible. Windy from the SW gusty at times',2,4],
['2015-12-15','Tuesday',40,46,10,16,14,3,5,26,47,'Fresh powder and groomed variable conditions','4WD or chains RECOMMENDED','16\" of fresh snow (19\" storm total; 31\" 72-hour total). NW wind 15 mph',2,3],
['2015-12-16','Wednesday',40,46,0,0,22,3,5,26,47,'Packed powder and groomed variable conditions','4WD or chains RECOMMENDED','Chilly with a high of 22 with N winds 6-10 mph',3,3],
['2015-12-17','Thursday',40,46,0,0,36,3,5,26,47,'Packed powder and groomed variable conditions','4WD or chains RECOMMENDED','high of 36 with N winds 5-11 mph breezy at times',3,3],
['2015-12-18','Friday',40,46,0,0,42,4,5,32,47,'Packed powder and groomed variable conditions','plowed and cindered','high of 42 with calm East winds 5-8 mph',3,1],
['2015-12-19','Saturday',40,46,0,0,44,5,5,35,47,'Packed powder and groomed variable conditions','Clear','high of 44 with SW winds 5-15 mph',3,0],
['2015-12-21','Monday',40,46,1,1,37,5,5,35,47,'Trace of New - Powder and packed powder','4WD or chains recommended','High of 37 NW winds 0-9 mph and sunny.',1,3],
['2015-12-22','Tuesday',40,46,0,0,34,5,5,35,47,'Packed powder variable conditions','4WD or chains recommended','high near 34. SW wind 17-22 mph. Chance of snow 70% with 1-2\" possible.',2,3],
['2015-12-23','Wednesday',38,40,1,2,32,5,5,35,47,'Packed powder variable conditions','4WD or chains required','high 32. with fog and W wind 10-17 mph. Chance of snow 50% with 1-2\" possible.',2,4],
['2015-12-24','Thursday',38,40,1,1,29,5,5,35,47,'Packed powder variable conditions','4WD or chains required','high 29. with 50% chance of snow wind 0-10 mph.',2,4],
['2015-12-25','Friday',38,40,1,1,18,5,5,35,47,'Packed powder variable conditions','4WD or chains required','high 18. Snow showers in the morning. 2-4\" possible. Windy gusty at times.',2,4],
['2015-12-26','Saturday',38,41,3,5,6,5,5,35,47,'Packed powder variable conditions','4WD or chains required','high of 6 NE wind 17-25 mph',2,4],
['2015-12-27','Sunday',38,41,0,0,28,5,5,35,47,'Packed powder variable conditions','4WD or chains recommended','high of 28 NE wind 8-15 mph Extended hrs 8:30-4:30',2,3],
['2015-12-28','Monday',38,41,0,0,24,5,5,35,47,'Packed powder variable conditions','plowed and cindered','high of 24 SW winds 5-15 mph',2,1],
['2015-12-29','Tuesday',38,41,1,1,18,5,5,35,47,'Packed powder variable conditions','Plowed and cindered','high of 18 W winds 5-9 mph',2,1],
['2015-12-30','Wednesday',38,41,0,0,24,5,5,35,47,'Packed powder variable conditions','Plowed and cindered','high of 24 W winds 5-12 mph',2,1],
['2015-12-31','Thursday',38,41,0,0,22,5,5,35,47,'Packed powder variable conditions','Clear','high of 22 E winds 0-7 mph',2,0],
['2016-01-01','Friday',38,41,0,0,27,5,5,35,47,'Packed powder variable conditions','Clear','high of 27 SE winds 9-11 mph',2,0],
['2016-01-02','Saturday',38,41,0,0,37,5,5,35,47,'Packed powder variable conditions','Clear','high of 37 S winds 8 mph',2,0],
['2016-01-03','Sunday',38,41,0,0,37,5,5,35,47,'Packed powder variable conditions','Clear','high of 37 S winds 5-12 mph',2,0],
['2016-01-04','Monday',38,41,0,0,32,5,5,35,47,'Packed powder variable conditions','Clear','high of 33 S winds 10-15 mph',2,0],
['2016-01-05','Tuesday',42,48,5,8,26,5,5,35,47,'Powder and Packed Powder conditions','4x4 or chains required','New Snow - 8\" with 2-4\" more today',1,4],
['2016-01-06','Wednesday',48,53,8,13,25,5,5,35,47,'Powder and Packed Powder conditions','4x4 or chains required','New Snow - 13\" with 2-4\" more today! Two day storm total so far is 21\"',1,4],
['2016-01-07','Thursday',60,64,12,18,17,4,5,35,47,'Powder and Packed Powder conditions','4x4 or chains required','New Snow - 18"" with 9-13\" more today! Two day storm total so far is 39\"',1,4],
['2016-01-08','Friday',68,72,10,16,17,4,5,35,47,'Powder and Packed Powder conditions','4x4 or chains required','New Snow - 16\" with 2-4\" more today! Four day storm total so far is at 55\"',1,4],
['2016-01-09','Saturday',70,70,1,3,24,5,5,35,47,'Powder and Packed Powder conditions','4x4 or chains recommended','New Snow 3\" Five day storm total 58\"',1,3],
['2016-01-10','Sunday',70,70,0,0,21,5,5,35,47,'Powder and Packed Powder conditions','Icy in some spots','1-3\" of snow expected today!  Five day storm total 58\"',1,2],
['2016-01-11','Monday',70,70,0,0,28,5,5,35,47,'Powder and Packed Powder conditions','Icy in some spots','Trace of new snow!',1,2],
['2016-01-12','Tuesday',70,70,0,0,33,5,5,35,47,'Powder and Packed Powder conditions','Icy in some spots','Blue Bird Day!',1,2],
['2016-01-13','Wednesday',70,70,0,0,40,5,5,35,47,'Powder and Packed Powder conditions','Icy in some spots','Blue Bird Day!',1,2],
['2016-01-14','Thursday',68,70,0,0,31,5,5,35,47,'Powder and Packed Powder conditions','Icy in some spots','Sunshine!  Awesome day!',1,2],
['2016-01-15','Friday',69,69,0,0,27,5,5,35,47,'Packed powder','Icy in some spots','3\"-5\" expected today!',2,2],
['2016-01-16','Saturday',68,68,4,6,28,5,5,35,47,'Packed powder','Icy in some spots','30% chance of snow today!',2,2],
['2016-01-17','Sunday',68,68,0,0,26,5,5,35,47,'Packed powder','Icy spots are present but no restrictions','Agassiz lift is on wind hold. All other lifts open until 4:30 p.m.',2,2],
['2016-01-18','Monday',68,68,0,0,31,5,5,35,47,'Packed powder','Icy spots are present but no restrictions','Sunny with a high near 43',2,2],
['2016-01-19','Tuesday',61,61,0,0,31,5,5,35,47,'Packed powder','Icy spots are present but no restrictions','Cloudy with a high near 38',2,2],
['2016-01-20','Wednesday',62,62,1,1,31,5,5,35,47,'Powder and Packed powder','Icy spots are present but no restrictions','Breezy with a high near 37',1,2],
['2016-01-22','Friday',62,62,0,0,44,5,5,39,47,'Packed powder','Icy spots are present but no restrictions','Mostly Sunny with a high near 44',2,2],
['2016-01-23','Saturday',62,62,0,0,32,5,5,39,47,'Packed powder','Icy spots are present but no restrictions','Mostly Sunny with Increasing Winds in the Afternoon',2,2],
['2016-01-25','Monday',62,62,0,0,24,5,5,36,47,'Packed powder','No restrictions','Mostly sunny!',2,0],
['2016-01-27','Wednesday',60,60,0,0,27,5,5,36,47,'Packed powder','No restrictions','Mostly sunny!',2,0],
['2016-01-28','Thursday',60,60,0,0,25,5,5,36,47,'Packed powder','No restrictions','Mostly sunny',2,0],
['2016-01-30','Saturday',60,60,0,0,24,5,5,36,47,'Packed powder','No restrictions','Wind Advisory in Effect from 10am to 7pm!',2,0],
['2016-02-01','Monday',69,69,10,13,20,5,5,38,47,'Powder and Packed powder','Chains and or 4x4 required','Heavy snow!',1,4],
['2016-02-03','Wednesday',70,70,0,0,11,5,5,38,47,'Powder and Packed Powder','Some Icy Spots','Sunny with moderate winds. Dress Warm!',1,2],
['2016-02-04','Thursday',70,70,0,0,12,5,5,38,47,'Powder and Packed Powder','Some Icy Spots','Sunny with moderate winds. Dress Warm!',1,2],
['2016-02-05','Friday',68,68,0,0,19,5,5,38,47,'Powder and Packed Powder','Some Icy Spots','Sunny and brisk with moderate winds. Dress Warm!',1,2],
['2016-02-07','Sunday',68,68,0,0,28,5,5,38,47,'Powder and Packed Powder','Some Icy Spots','Sunny with light winds',1,2],
['2016-02-12','Friday',60,60,0,0,38,5,5,38,47,'Spring like conditions','Clear','Warm and sunny',5,0],
['2016-02-14','Sunday',58,58,0,0,38,5,5,38,47,'Spring like conditions','Clear','Warm and sunny moderate winds in the afternoon.',5,0],
['2016-02-15','Monday',58,58,0,0,28,5,5,38,47,'Spring like conditions','Clear','Blustery day winds 5-10mph decreasing by afternoon.',5,0],
['2016-02-16','Tuesday',56,56,0,0,34,5,5,38,47,'Spring like conditions','Clear','Sunny Winds Increasing By Afternoon.',5,0],
['2016-02-17','Wednesday',56,56,0,0,34,5,5,38,47,'Spring like conditions','Clear','Sunny Winds Increasing By Afternoon.',5,0],
['2016-02-19','Friday',57,57,0,0,30,5,5,38,47,'Spring like conditions','Clear','Mostly sunny and breezy',5,0],
['2016-02-20','Saturday',57,57,0,0,30,5,5,38,47,'Variable Spring-like conditions','Clear','Mostly Sunny',5,0],
['2016-02-21','Sunday',57,57,0,0,34,5,5,38,47,'Variable Spring-like conditions','Clear','Beautiful and Sunny Light winds.',5,0],
['2016-02-24','Wednesday',56,56,0,0,38,5,5,38,47,'Variable Spring-like conditions','Clear','Sunny with a slight breeze before noon',5,0],
['2016-02-25','Thursday',55,55,0,0,37,5,5,38,47,'Variable Spring-like conditions','Clear','Sunny with a slight breeze!',5,0],
['2016-02-26','Friday',57,57,0,0,43,5,5,38,47,'Variable Spring-like conditions','Clear','Sunny!',5,0],
['2016-02-27','Saturday',57,57,0,0,37,5,5,38,47,'Variable Spring-like conditions','Clear','Sunny!',5,0],
['2016-02-28','Sunday',57,57,0,0,40,5,5,39,47,'Variable Spring-like conditions','Clear','Sunny!',5,0],
['2016-02-29','Monday',57,57,0,0,36,4,5,39,47,'Variable Spring-like conditions','Clear','Mostly Sunny with a Chance of Light Snow After 11am',5,0],
['2016-03-01','Tuesday',57,57,0,0,36,4,5,39,47,'Variable Spring-like conditions','Clear','Mostly Sunny with a Chance of Light Snow After 11am',5,0],
['2016-03-02','Wednesday',52,52,0,0,32,4,5,38,47,'Variable Spring-like conditions','Clear','Mostly Sunny with light winds 15-20 mph',5,0],
['2016-03-03','Thursday',52,52,0,0,40,4,5,36,47,'Variable Spring-like conditions','Clear','Mostly Sunny with light winds 15-20 mph',5,0],
['2016-03-04','Friday',52,52,0,0,35,4,5,37,47,'Variable Spring-like conditions','Clear','Mostly Sunny with winds 10-20 mph',5,0],
['2016-03-06','Sunday',52,52,0,0,33,4,5,37,47,'Variable Spring-like conditions','Clear','Light snow throughout the day with winds10-20 mph from the WSW',5,0],
['2016-03-07','Monday',53,53,3,3,30,4,5,37,47,'Powder packed powder','chains and or 4 wheel drive is recommended','Light snow throughout the day with winds15-20 mph from the WSW',2,3],
['2016-03-08','Tuesday',56,56,1,1,30,4,5,37,47,'Powder packed powder','chains and or 4 wheel drive is recommended','Light snow throughout the day with winds15-25 mph from the NE',2,3],
['2016-03-09','Wednesday',54,54,0,0,24,4,5,37,47,'Powder packed powder','Clear','Sunny and breezy. Winds15-25 mph from the north northeast with gusts as high as 31 mph',2,0],
['2016-03-10','Thursday',51,51,0,0,31,4,5,37,47,'Powder packed powder','Clear','Sunny and breezy. Slight winds 0-5 mph from the north northeast.',2,0],
['2016-03-11','Friday',51,51,0,0,33,4,5,37,47,'Powder packed powder','Clear','Sunny and breezy! Winds 15-25 mph from the south southwest gusting up to 35 mph.',2,0],
['2016-03-13','Sunday',55,55,3,3,24,4,5,37,47,'Powder packed powder','Plowed CHAINS or 4WD Recommended','Light snow in the morning until 11 a.m.! Sunny by afternoon and breezy!',2,3],
['2016-03-14','Monday',54,54,0,0,30,4,5,36,47,'Powder packed powder','Clear','Sunny and breezy. Winds out of the west at 15-25 mph with gusts up to 30 mph.',2,0],
['2016-03-15','Tuesday',54,54,0,0,30,4,5,36,47,'Powder packed powder','Clear','Sunny and breezy. Winds out of the NW 10-20 mph.',2,0],
['2016-03-16','Wednesday',54,54,0,0,32,4,5,36,47,'Powder packed powder','Clear','Sunny. Light winds out of the north northeast around 6 mph.',2,0],
['2016-03-17','Thursday',52,52,0,0,32,4,5,36,47,'Spring Conditions','Clear','Sunny. Light winds out of the north northeast around 6 mph.',5,0],
['2016-03-18','Friday',52,52,0,0,42,4,5,36,47,'Spring Conditions','Clear',"Sunny. Light winds out of the north northeast around 6 mph. It's going to be a great day!",5,0],
['2016-03-19','Saturday',52,52,0,0,41,4,5,36,47,'Spring Conditions','Clear',"Sunny. Light variable winds. It's going to be a great day!",5,0],
['2016-03-20','Sunday',52,52,0,0,41,4,5,36,47,'Spring Conditions','Clear',"Sunny. Light variable winds. It's going to be a great day!",5,0],
['2016-03-21','Monday',52,52,0,0,44,3,5,36,47,'Spring Conditions','Clear','Sunny. Breezy with winds out of the southwest at 15-30 mph with gusts up to 40 mph.',5,0],
['2016-03-28','Monday',50,50,0,0,30,0,5,32,47,'Spring Conditions','Clear','Sunny and Breezy! Chance of snow for the remainder of the week!!!!!',5,0],
['2016-03-29','Tuesday',50,50,0,0,30,0,5,32,47,'Spring Conditions','Clear','Sunny and Breezy! Chance of snow for the remainder of the week!!!!!',5,0]]

#print(len(olddata),len(olddata[0]),len(oldcols),len(data),len(data[0]),len(cols))

oldDate = [0 for i in range(len(olddata))]
oldMinDepth = [0 for i in range(len(olddata))]
oldMaxDepth = [0 for i in range(len(olddata))]
oldNewBaseSnow = [0 for i in range(len(olddata))]
oldNewPeakSnow = [0 for i in range(len(olddata))]
oldSnowHazard = [0 for i in range(len(olddata))]
oldRoadHazard = [0 for i in range(len(olddata))]
oldOpenLifts = [0 for i in range(len(olddata))]
oldOpenTrails = [0 for i in range(len(olddata))]
oldTemp = [0 for i in range(len(olddata))]

Date = [0 for i in range(len(data))]
MinDepth = [0 for i in range(len(data))]
MaxDepth = [0 for i in range(len(data))]
NewBaseSnow = [0 for i in range(len(data))]
NewPeakSnow = [0 for i in range(len(data))]
SnowHazard = [0 for i in range(len(data))]
RoadHazard = [0 for i in range(len(data))]
OpenLifts = [0 for i in range(len(data))]
OpenTrails = [0 for i in range(len(data))]
AzSbTemp = [0 for i in range(len(data))]
LowTemp = [0 for i in range(len(data))]
HighTemp = [0 for i in range(len(data))]

def date_to_num(datestr):
    #datestr must be in format YYYY-MM-DD
    #returns format YYMMDD
    #datestr.split(sep='-') doesn't work
    #return int(datestr[0:4]+datestr[5:7]+datestr[8:]) #gets formatted as sci-notation
    #return int(datestr[2:4]+datestr[5:7]+datestr[8:])  #jumps at turn of year
    #This is strictly for snow season purposes
    yr = int(datestr[0:4])
    mo = int(datestr[5:7])
    dy = int(datestr[8:])
    
    ndays = {1:31,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
    #its a leapyear
    if yr%4 != 0:
        ndays[2]=28
    elif yr%100 != 0:
        ndays[2]=29
    elif yr%400 != 0:
        ndays[2]=28
    else:
        ndays[2]=29
    
    datenum = (mo - 1) + dy/(ndays[mo]+1)
        
    if mo > 7:
        datenum -= 12
        
    return datenum

def wrangle_data():
    global oldDate, oldMinDepth, oldMaxDepth,oldNewBaseSnow,oldNewPeakSnow
    global oldSnowHazard,oldRoadHazard,oldOpenLifts,oldOpenTrails,oldTemp
    global Date,MinDepth,MaxDepth,NewBaseSnow,NewPeakSnow,SnowHazard
    global RoadHazard,OpenLifts,OpenTrails,AzSbTemp,LowTemp,HighTemp
    #makes data easier to use for glowscript plotting
    #subplot: 0,0:2, (all slice ends exclusive)
    #'MinDepth','MaxDepth','NewBaseSnow','NewPeakSnow'
    #subplot: 1,0:2
    #'SnowHazard','RoadHazard'
    #subplot: 2,0:2
    #'OpenLifts','OpenTrails'
    #subplot: 3,0:3
    #olddata: 'Temp'
    #newdata: 'AzSbTemp','LowTemp','HighTemp'
    #ALL olddata has prefix old: e.g. oldMinDepth, oldTemp
    #olddate vs date
    
    for i in range(len(olddata)):
        oldDate[i] = olddata[i][0]
        oldMinDepth[i] = [date_to_num(olddata[i][0]),olddata[i][2]]
        oldMaxDepth[i] = [date_to_num(olddata[i][0]),olddata[i][3]]
        oldNewBaseSnow[i] = [date_to_num(olddata[i][0]),olddata[i][4]]
        oldNewPeakSnow[i] = [date_to_num(olddata[i][0]),olddata[i][5]]
        oldSnowHazard[i] = [date_to_num(olddata[i][0]),olddata[i][-2]]
        oldRoadHazard[i] = [date_to_num(olddata[i][0]),olddata[i][-1]]
        oldOpenLifts[i] = [date_to_num(olddata[i][0]),olddata[i][7]]
        oldOpenTrails[i] = [date_to_num(olddata[i][0]),olddata[i][9]]
        oldTemp[i] = [date_to_num(olddata[i][0]),olddata[i][6]]
        
    for i in range(len(data)):
        Date[i] = data[i][0]
        MinDepth[i] = [date_to_num(data[i][0]),data[i][2]]
        MaxDepth[i] = [date_to_num(data[i][0]),data[i][3]]
        NewBaseSnow[i] = [date_to_num(data[i][0]),data[i][4]]
        NewPeakSnow[i] = [date_to_num(data[i][0]),data[i][5]]
        SnowHazard[i] = [date_to_num(data[i][0]),data[i][-2]]
        RoadHazard[i] = [date_to_num(data[i][0]),data[i][-1]]
        OpenLifts[i] = [date_to_num(data[i][0]),data[i][10]]
        OpenTrails[i] = [date_to_num(data[i][0]),data[i][12]]
        AzSbTemp[i] = [date_to_num(data[i][0]),data[i][6]]
        LowTemp[i] = [date_to_num(data[i][0]),data[i][7]]
        HighTemp[i] = [date_to_num(data[i][0]),data[i][8]]
    

wrangle_data()

cw = 800
ch = cw

gw = 0.7*cw
gh = 0.333*gw

#subplot: 0,0:2, (all slice ends exclusive)
#'MinDepth','MaxDepth','NewBaseSnow','NewPeakSnow'
#subplot: 1,0:2
#'SnowHazard','RoadHazard'
#subplot: 2,0:2
#'OpenLifts','OpenTrails'
#subplot: 3,0:3
#olddata: 'Temp'
#newdata: 'AzSbTemp','LowTemp','HighTemp'
#ALL olddata has prefix old: e.g. oldMinDepth, oldTemp
#olddate vs date
clr = [color.blue,color.green,color.red,color.black,color.magenta,color.orange,color.yellow]
title1='Snowbowl stats 2015-16 (left), 2016-17 (right); x=0 at 01 Jan'

gInch00 = graph(width=gw,height=gh,align='left',ymin=0,title=title1)
gInch01 = graph(width=gw,height=gh,align='right',ymin=0)
gHzd10 = graph(width=gw,height=gh,align='left',ymin=0,ymax=6)
gHzd11 = graph(width=gw,height=gh,align='right',ymin=0,ymax=6)
gOpen20 = graph(width=gw,height=gh,align='left',ymin=0,ymax=50)
gOpen21 = graph(width=gw,height=gh,align='right',ymin=0,ymax=50)
gTemp30 = graph(width=gw,height=gh,align='left')
gTemp31 = graph(width=gw,height=gh,align='right')

#plot them
p00a = gcurve(graph=gInch00,pos=oldMinDepth,label='MinDepth (in)',color=clr[0])
p00b = gcurve(graph=gInch00,pos=oldMaxDepth,label='MaxDepth (in)',color=clr[1])
p00c = gcurve(graph=gInch00,pos=oldNewBaseSnow,label='NewBaseSnow (in)',color=clr[2])
p00d = gcurve(graph=gInch00,pos=oldNewPeakSnow,label='NewPeakSnow (in)',color=clr[3])
#p00a.plot(oldMinDepth)
#p00b.plot(oldMaxDepth)
#p00c.plot(oldNewBaseSnow)
#p00d.plot(oldNewPeakSnow)
#
p01a = gcurve(graph=gInch01,pos=MinDepth,label='MinDepth (in)',color=clr[0])
p01b = gcurve(graph=gInch01,pos=MaxDepth,label='MaxDepth (in)',color=clr[1])
p01c = gcurve(graph=gInch01,pos=NewBaseSnow,label='NewBaseSnow (in)',color=clr[2])
p01d = gcurve(graph=gInch01,pos=NewPeakSnow,label='NewPeakSnow (in)',color=clr[3])
#
p10a = gcurve(graph=gHzd10,pos=oldSnowHazard,label='SnowHazard [0-5]',color=clr[3])
p10b = gcurve(graph=gHzd10,pos=oldRoadHazard,label='RoadHazard [0-4]',color=clr[2])
#
p11a = gcurve(graph=gHzd11,pos=SnowHazard,label='SnowHazard [0-5]',color=clr[3])
p11b = gcurve(graph=gHzd11,pos=RoadHazard,label='RoadHazard [0-4]',color=clr[2])
#
p20a = gcurve(graph=gOpen20,pos=oldOpenLifts,label='OpenLifts',color=clr[0])
p20b = gcurve(graph=gOpen20,pos=oldOpenTrails,label='OpenTrails',color=clr[1])
#
p21a = gcurve(graph=gOpen21,pos=OpenLifts,label='OpenLifts',color=clr[0])
p21b = gcurve(graph=gOpen21,pos=OpenTrails,label='OpenTrails',color=clr[1])
#
p30a = gcurve(graph=gTemp30,pos=oldTemp,label='snowbowlTemp (F)',color=clr[2])
#
p31a = gcurve(graph=gTemp31,pos=AzSbTemp,label='snowbowlTemp (F)',color=clr[2])
p31b = gcurve(graph=gTemp31,pos=LowTemp,label='LowTemp (F)',color=clr[0])
p31c = gcurve(graph=gTemp31,pos=HighTemp,label='HighTemp (F)',color=clr[1])


