monthToDays = {
        1: 0,
        2: 31,
        3: 60,
        4: 91,
        5: 121,
        6: 152,
        7: 182,
        8: 213,
        9: 244,
        10: 274,
        11: 305,
        12: 335
}





startTime = []
endTime = []
startHours = []
endHours = []
startLat = []
startLon = []
endLat = []
endLon = []
lines = 0

outfile = open("taxi.txt","w")
for line in open("train.csv","r"):
    tokens = line.split(",")
    startTime.append(tokens[2])
    endTime.append(tokens[3])
    startLat.append(tokens[6])
    startLon.append(tokens[5])
    endLat.append(tokens[8])
    endLon.append(tokens[7])
    lines += 1

for point in startTime:
    hours = 0
    date, time = point.split(" ")
    year, month, day = date.split("-")
    hour, minute, second = time.split(":")


    hours += monthToDays[float(month)] * 24
    hours += float(day) * 24
    hours += float(hour)
    hours += float(minute)/60
    hours += float(second)/(3600)
    startHours.append(hours)

for point in endTime:
    hours = 0
    date, time = point.split(" ")
    year, month, day = date.split("-")
    hour, minute, second = time.split(":")

    hours += monthToDays[float(month)] * 24
    hours += float(day) * 24
    hours += float(hour)
    hours += float(minute)/60
    hours += float(second)/(3600)
    endHours.append(hours)
    
for i in range(lines):
    assert float(startHours[i]) < float(endHours[i])
    outfile.write(str(i) + " 1 1 " + str(startHours[i]) + " " + str(endHours[i]) + " " + str(endHours[i]-startHours[i]) + " " + startLat[i] + " " + startLon[i] + " " + endLat[i] + " " + endLon[i] + "\n")

outfile.close()

