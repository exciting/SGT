from sgt import  *
distortions=loadObject("distwenergy")
#distortions.make_report()

#distortions.select_data()

elasticcontants=distortions.get_elastic_c(selectdata=[0.02,2] )
print elasticcontants
distortions.make_report()
