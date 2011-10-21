from stool import  *
distortions=loadObject("distwenergy")
distortions.make_report()
distortions.select_data()

elasticcontants=distortions.get_elastic_c()

distortions.make_report()