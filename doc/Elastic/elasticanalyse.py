from sgt import  *
distortions=loadObject("distwenergy")
#distortions.make_report()

distortions.show_data()

elasticcontants=distortions.get_elastic_c(selectdata=[0.02,2] )
distortions.make_report()
