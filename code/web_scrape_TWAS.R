
#########################  Web scraping TWAS hub  #############################

library('rvest')
webpage <- html_session('http://twas-hub.org/traits/')
webpage %>% html_nodes("div.main > a") %>% html_attr('href')


my_session <- html_session("https://scrapethissite.com/pages/simple/")
my_nodes <- my_session %>% html_nodes(".country")
my_texts <- my_session %>% html_nodes(".country-capital") %>% html_text()


library(splashr)
#> Warning: le package 'splashr' a été compilé avec la version R 3.4.4
sp <- splash("192.168.99.100")
page <- render_html(sp, url = "https://www.nhl.com/gamecenter/phi-vs-bos/1974/05/07/1973030311#game=1973030311,game_state=final")

library(rvest)
#> Le chargement a nécessité le package : xml2
page %>%
  html_nodes("div.name > strong > a") %>%
  html_attr("href")



library(dplyr)
library(rvest)

url<-html("beer.html")
selector_name<-".brewery"
fnames<-html_nodes(x = url, css = selector_name) %>%
  html_text()
head(fnames)
fnames


url <- paste0(
  "https://web.archive.org/web/20190202054736/",
  "https://www.boxofficemojo.com/movies/?id=ateam.htm"
)
ateam <- read_html(url)
html_nodes(ateam, "center")
html_nodes(ateam, "center font")
html_nodes(ateam, "center font b")


library('rvest')
webpage <- html_session('http://twas-hub.org/traits/')
ateam <- read_html(webpage)
tt = html_nodes(ateam, "body a") %>% html_attr('href')
traits = tt[grep("/traits", tt)][-1]
data = tt[grep("/data", tt)]
description = html_nodes(ateam, "body a") %>% html_text()
description = description[-(1:6)]
descr2 = description[c(TRUE, FALSE)]
data = paste0(" http://twas-hub.org", data)
traits = paste0(" http://twas-hub.org/", traits)

final_df = cbind.data.frame(traits, data, descr2)
colnames(final_df) = c("Trait", "Data_path", "Description")
write.table(final_df, file = "/n/groups/price/kushal/singlecellLDSC/data/TWAS-Hub/TWAS_models_traits.txt",
            row.names = F, quote=F, sep = "\t")



