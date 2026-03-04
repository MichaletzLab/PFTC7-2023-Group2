at.subset3 <- read.csv("data/raw.discardHooks_data.csv")
at.subset3 <- at.subset3 %>%
  mutate(iWUE = A/gsw)%>%
  mutate(eWUE = A/E)

species.key = read.csv("data/Faster_Key.csv")%>%
  rename(curveID=Obs)%>%
  rename(site=SiteID)
species.key.n=read.csv("data/Norway.Key.csv")%>%
  mutate(site=site+5)%>%
  rename(Elevation = Elevation.masl)

Sfield.results_meta_discard.hooks <- read.csv("schoolfield.discard.hooks.SANW.csv")
gsw.Sfield.results_meta_discard.hooks <- read.csv("schoolfield.gsw.discard.hooks.SANW.csv")
iWUE.Sfield.results_meta_discard.hooks <- read.csv("schoolfield.iWUE.discard.hooks.SANW.csv")
eWUE.Sfield.results_meta_discard.hooks <- read.csv("schoolfield.eWUE.discard.hooks.SANW.csv")
