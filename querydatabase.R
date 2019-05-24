#Automatically query the database

#I need to use 32-bit R if I want to connect to the database automatically. 

#Set up a function to get data from the FRP database
#"path" should be the location and name of the database
#"type" should be what kind of data you want. Either "zoop", "fish", "inverts", "phyto", or "WQ", or "stations"

GetFRPdata = function(path, type) {
  library(RODBC)
  channel <- odbcConnectAccess2007(path)
  
  if (type == "inverts") {
 data =   sqlQuery(channel, query = paste("SELECT SampleInfo.SampleName AS SampleID, 
                                    SiteVisit.Date, SampleInfo.Starttime AS [time], SampleInfo.[Sample type] AS Sampletype, 
                                    SampleInfo.VolumeEstimated AS Volume, SiteVisit.Station, SampleInfo.LatStart AS Latitude, 
                                    [LongStart]*-1 AS Longitude, SiteVisit.Temp, SiteVisit.SC, SiteVisit.pH, SiteVisit.DO, 
                                    SiteVisit.Secchi, SiteVisit.Turbidity, SiteVisit.Chlorophyll, SiteVisit.FDOM, 
                                    SiteVisit.Tide, SiteVisit.Microcystis, Inverts.subsampled, SpeciesCodes.CommonName, 
                                    [Invert Catch].presence, [Invert Catch].TotalCount, 
                                    [Invert Catch].[TotalCount]/([Inverts].[subsampled]/100) AS AdjCount
                                    FROM SpeciesCodes INNER JOIN 
                                    (((SiteVisit INNER JOIN SampleInfo ON SiteVisit.VisitNo = SampleInfo.[Site visit]) 
                                    INNER JOIN Inverts ON SampleInfo.SampleName = Inverts.[miSample ID]) INNER JOIN [Invert Catch] 
                                    ON (SampleInfo.SampleName = [Invert Catch].SampleName) AND (Inverts.[miSample ID] = 
                                    [Invert Catch].SampleName)) ON SpeciesCodes.[Zoo Code] = [Invert Catch].InvertCode
                                    WHERE (((SampleInfo.[Sample type])='Mysid net (no sled)' Or 
                                    (SampleInfo.[Sample type])='larval sled benthic' Or 
                                    (SampleInfo.[Sample type])='larval sled oblique' Or 
                                    (SampleInfo.[Sample type])='neuston trawl' Or 
                                    (SampleInfo.[Sample type])='sweep net' Or 
                                    (SampleInfo.[Sample type])='PCV core' Or 
                                    (SampleInfo.[Sample type])='Ponar grab' Or 
                                    (SampleInfo.[Sample type])='Petite Ponar'));"))
 return(data)
  }
  
    else if (type == "zoop") {
      data = sqlQuery(channel, query = paste("SELECT SampleInfo.SampleName AS SampleID, SiteVisit.Date, 
SampleInfo.Starttime AS [time], SampleInfo.[Sample type] AS Sampletype, 
SampleInfo.VolumeEstimated AS Volume, SiteVisit.Station, SampleInfo.LatStart AS Latitude, 
[LongStart]*-1 AS Longitude, SiteVisit.Temp, SiteVisit.SC, SiteVisit.pH, SiteVisit.DO, 
SiteVisit.Secchi, SiteVisit.Turbidity, SiteVisit.Chlorophyll, SiteVisit.FDOM, SiteVisit.Tide, 
SiteVisit.Microcystis, Zooplankton.Dilution, Zooplankton.CellNumber, SpeciesCodes.CommonName, [Zoo Catch].Count
                                             FROM ((SiteVisit INNER JOIN SampleInfo ON 
SiteVisit.VisitNo = SampleInfo.[Site visit]) INNER JOIN Zooplankton ON 
SampleInfo.SampleName = Zooplankton.[zSample ID]) INNER JOIN (SpeciesCodes INNER JOIN [Zoo Catch] 
ON SpeciesCodes.[Zoo Code] = [Zoo Catch].ZooCode) ON Zooplankton.ZoopID = [Zoo Catch].ZoopID
                                             WHERE (((SampleInfo.[Sample type])='Zoopnet (no sled)' 
Or (SampleInfo.[Sample type])='Zoopnet oblique sled' Or 
(SampleInfo.[Sample type])='Zoopnet benthic sled'));"))
 return(data)
      }
else if (type == "stations"){
  data = sqlFetch(channel, "Stations2")
  return(data)
}
  
  else if (type == "phytoplankton"){
    data =   sqlQuery(channel, query = paste( "SELECT SiteVisit.Station, Stations2.Region, 
SiteVisit.Date, SiteVisit.Temp, SiteVisit.SC, SiteVisit.pH, SiteVisit.DO, SiteVisit.Secchi, 
SiteVisit.Turbidity, SiteVisit.Chlorophyll, SiteVisit.FDOM, SiteVisit.PC, SiteVisit.Microcystis, 
SampleInfo.[Habitat type], SampleInfo.SampleName, phytoplankton.Volume, phytoplankton.PercentCounted, 
phytoplankton.Taxon, phytoplankton.CellCount, phytoplankton.CellspermL, SampleInfo.Comments,  phytoplankton.BiovolumeperuL
FROM Stations2 INNER JOIN ((SiteVisit INNER JOIN SampleInfo ON SiteVisit.VisitNo = SampleInfo.[Site visit]) 
                                              INNER JOIN phytoplankton ON SampleInfo.SampleName = phytoplankton.SampleID) ON Stations2.Station = SiteVisit.Station;"))
    return(data)
  }
  
  
 
  else if (type == "chlorophyll"){
    data =   sqlQuery(channel, query = paste( "SELECT SiteVisit.Station, SiteVisit.Date, SiteVisit.Chlorophyll, 
Stations2.Region, SiteVisit.VisitNo
FROM Stations2 INNER JOIN SiteVisit ON Stations2.Station = SiteVisit.Station
                                              WHERE (((SiteVisit.Chlorophyll) Is Not Null));"))
    return(data)
  }
  
  
  else print("Sorry, haven't gotten that far yet")
  
  odbcClose(channel)
                                             
  
}




