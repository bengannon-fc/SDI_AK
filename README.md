# Suppression Difficulty Index (SDI) for Alaska

Suppression Difficulty Index (SDI) is a spatial decision support product developed by a group of wildland firefighting experts in Spain (Rodríguez y Silva et al. 2014, 2020) that has been adapted for use in the United States by the USDA Forest Service Rocky Mountain Research Station. The formulation for SDI is premised on the concept that fire behavior acts as a driving force in the numerator and conditions that facilitate suppression act as resisting forces in the denominator (Equation 1). The numerator, also called the energy behavior index (_Ice_), grows in response to conditions that increase fire behavior, so SDI will increase where and when fuel, topography, and weather conditions promote extreme fire behaviors. The denominator sums the values for sub-indices representing how accessibility (_Ia_), penetrability (_Ip_), and fireline creation (_Ic_) influence suppression opportunities across the landscape. Higher values for sub-indices in the denominator indicate better conditions for suppression, which lowers SDI for a given energy behavior index value in the numerator.   

**Equation 1**  _SDI = Ice/(Ia + Ip + Ic + 1)_

Feedback from early users in the United States suggested that SDI underestimated the difficulty of suppressing fire on steep slopes, especially ones with sparse fuels that were mapped as non-burnable by LANDFIRE. While steep slopes with little fuel (e.g., scree slopes, cliff bands, etc.) may act as natural barriers to fire, these areas are not ideal locations to deploy hand crews or machines to construct containment line. Starting in the 2021 season, an additional slope factor was applied to increase SDI on steep slopes regardless of fire behavior or other operational factors.

**Equation 2**  _SDI = Ice/(Ia + Ip + Ic + 1) + Slope Factor_

SDI provides a snapshot of how potential fire behavior and operational factors combine to influence suppression difficulty across broad landscapes. SDI tends to align well with firefighter expectations for low suppression difficulty near major roads, significant waterbodies, and flatter terrain, and for high suppression difficulty in steep, forested areas without recent history of wildfire or fuel treatment. SDI is moderately sensitive to the fire behavior scenario used to model the energy behavior index (numerator); SDI will rise with increasingly extreme fire weather inputs, but the spatial patterns and their interpretations for fire management do not drastically change across the range of typical large fire weather conditions. Hence, it is usually not necessary to model SDI for incident specific fire weather conditions.

Full details on the standard SDI calculation methods can be found in Risk Management Assistance (RMA) SharePoint site.

**Alaska Modifications**

The Alaska Wildland Fire Coordinating Group's Fire Modeling and Analysis Committee adopted the following modifications to the standard SDI calculation methods:

1) The standard SDI calculations include an influence of slope aspect on the penetrability sub index to reflect that sun exposed south facing slopes are harder to work on than shaded north facing slopes ignoring vegetation and fuel differences. Aspect is not a strong influence on firefighting difficulty in Alaska except through its influence on vegetation and fuels, which are already captured in other model components. The aspect component was removed from the model and the pentrability sub-index denominator was changed from a value of 5 to a value of 4 to reflect one less factor in the sub-index numerator.

2)  
