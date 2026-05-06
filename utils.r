# Creating a mapping from ISO3 codes to IHME regions

# Create the dataframe with ISO3 codes and corresponding IHME regions
iso3_ihme_mapping <- data.frame(
  iso3 = c(
    # East Asia
    "HKG", "CHN", "TWN", "PRK", "JPN", "KOR", "MNG",
    "KOR",
    "KOR",
    
    # Southeast Asia
    "KHM", "IDN", "LAO", "MYS", "MMR", "PHL", "SGP", "THA", "TLS", "VNM", "BRN", "MDV",
    
    # Oceania
    "ASM", "AUS", "FJI", "FSM", "GUM", "KIR", "MHL", "NRU", "NCL", "NZL", "NIU", 
    "MNP", "PLW", "PNG", "WSM", "SLB", "TON", "TUV", "VUT", "COK", "PYF",
    
    # Central Asia
    "KAZ", "KGZ", "TJK", "TKM", "UZB", "GEO", "ARM", "AZE",
    "KGZ",  # Kyrgyzstan is also known as Kyrgyz Republic
    
    # Central Europe
    "ALB", "BGR", "BIH", "HRV", "CZE", "HUN", "MKD", "MNE", "POL", "ROU", "SRB", "SVK", "SVN",
    "CZE",  # Czech Republic is also known as Czechia
    "SVK",  # Slovakia is also known as Slovak Republic
    "MKD",  # North Macedonia is also known as Macedonia
    
    # Eastern Europe
    "BLR", "EST", "LVA", "LTU", "MDA", "RUS", "UKR",
    
    # Western Europe
    "AND", "AUT", "BEL", "CYP", "DNK", "FIN", "FRA", "DEU", "GRC", "ISL", "IRL", "ITA", 
    "LUX", "MLT", "MCO", "NLD", "NOR", "PRT", "SMR", "ESP", "SWE", "CHE", "GBR",
    
    # High-income North America
    "CAN", "USA",
    
    # Caribbean
    "ATG", "BHS", "BRB", "CUB", "DOM", "GRD", "HTI", "JAM", "PRI", "LCA", "VCT", "TTO", "VGB", "VIR", "DMA", "KNA", "GUY", "SUR",
    "LCA",  # Saint Lucia is also known as St. Lucia
    "VCT",  # Saint Vincent and the Grenadines is also known as St. Vincent and the Grenadines
    
    # Andean Latin America
    "BOL", "ECU", "PER",
    
    # Central Latin America
    "COL", "CRI", "SLV", "GTM", "HND", "MEX", "NIC", "PAN", "VEN", "BLZ",
    
    # Southern Latin America
    "ARG", "CHL", "URY",
    
    # Tropical Latin America
    "BRA", "PRY",
    
    # North Africa and Middle East
    "AFG", "DZA", "BHR", "EGY", "IRN", "IRQ", "JOR", "KWT", "LBN", "LBY", "MAR", "OMN", 
    "PSE", "QAT", "SAU", "SDN", "SYR", "TUN", "TUR", "ARE", "YEM", "ISR",
    
    # South Asia
    "BGD", "BTN", "IND", "NPL", "PAK", "LKA",
    
    # Central Sub-Saharan Africa
    "AGO", "CAF", "COG", "COD", "GNQ", "GAB",
    "COD",  # Democratic Republic of the Congo is also known as Congo, Dem. Rep.
    "COG",  # Republic of the Congo is also known as Congo, Rep.
    
    # Eastern Sub-Saharan Africa
    "BDI", "COM", "DJI", "ERI", "ETH", "KEN", "MDG", "MWI", "MOZ", "RWA", 
    "SOM", "SSD", "TZA", "UGA", "ZMB", "MUS", "SYC",
    
    # Southern Sub-Saharan Africa
    "BWA", "LSO", "NAM", "ZAF", "SWZ", "ZWE", "ESH",
    
    # Western Sub-Saharan Africa
    "BEN", "BFA", "CMR", "CPV", "TCD", "CIV", "GMB", "GHA", "GIN", "GNB", 
    "LBR", "MLI", "MRT", "NER", "NGA", "SEN", "SLE", "TGO", "STP",
    "CIV",  # Ivory Coast is also called "Cote d'Ivoire" in some contexts
    "CIV",  # Ivory Coast is also called "Cote d'Ivoire" in some contexts
    "GMB",   # Gambia is also referred to as "Gambia, The"
    "STP"   # Sao Tome and Principe without accents
  ),
  country_name = c(
    # East Asia
    "Hong Kong","China", "Taiwan", "Korea, North", "Japan", "Korea, South", "Mongolia",
    "Korea, Rep.", 
    "South Korea",
    # Southeast Asia
    "Cambodia", "Indonesia", "Laos", "Malaysia", "Myanmar", "Philippines", "Singapore", "Thailand", "Timor-Leste", "Vietnam", "Brunei", "Maldives",
    # Oceania
    "American Samoa", "Australia", "Fiji", "Micronesia", "Guam", "Kiribati", "Marshall Islands", "Nauru", "New Caledonia", "New Zealand", "Niue",
    "Northern Mariana Islands", "Palau", "Papua New Guinea", "Samoa", "Solomon Islands", "Tonga", "Tuvalu", "Vanuatu", "Cook Islands", "French Polynesia",
    # Central Asia
    "Kazakhstan", "Kyrgyzstan", "Tajikistan", "Turkmenistan", "Uzbekistan", "Georgia", "Armenia", "Azerbaijan",
    "Kyrgyz Republic",
    # Central Europe
    "Albania", "Bulgaria", "Bosnia and Herzegovina", "Croatia", "Czech Republic", "Hungary", "North Macedonia", "Montenegro", "Poland", "Romania", "Serbia", "Slovakia", "Slovenia",
    "Czechia",
    "Slovak Republic",
    "Macedonia",
    # Eastern Europe
    "Belarus", "Estonia", "Latvia", "Lithuania", "Moldova", "Russia", "Ukraine",
    # Western Europe
    "Andorra", "Austria", "Belgium", "Cyprus", "Denmark", "Finland", "France", "Germany", "Greece", "Iceland", "Ireland", "Italy",
    "Luxembourg", "Malta", "Monaco", "Netherlands", "Norway", "Portugal", "San Marino", "Spain", "Sweden", "Switzerland", "United Kingdom",
    # High-income North America
    "Canada", "United States",
    # Caribbean
    "Antigua and Barbuda", "Bahamas", "Barbados", "Cuba", "Dominican Republic", "Grenada", "Haiti", "Jamaica", "Puerto Rico", "Saint Lucia", "Saint Vincent and the Grenadines", "Trinidad and Tobago", "British Virgin Islands", "US Virgin Islands", "Dominica", "Saint Kitts and Nevis", "Guyana", "Suriname",
    "St. Lucia",
    "St. Vincent and the Grenadines",
    # Andean Latin America
    "Bolivia", "Ecuador", "Peru",
    # Central Latin America
    "Colombia", "Costa Rica", "El Salvador", "Guatemala", "Honduras", "Mexico", "Nicaragua", "Panama", "Venezuela", "Belize",
    # Southern Latin America
    "Argentina", "Chile", "Uruguay",
    # Tropical Latin America
    "Brazil", "Paraguay",
    # North Africa and Middle East
    "Afghanistan", "Algeria", "Bahrain", "Egypt", "Iran", "Iraq", "Jordan", "Kuwait", "Lebanon", "Libya", "Morocco", "Oman",
    "Palestine", "Qatar", "Saudi Arabia", "Sudan", "Syria", "Tunisia", "Turkey", "United Arab Emirates", "Yemen", "Israel",
    # South Asia
    "Bangladesh", "Bhutan", "India", "Nepal", "Pakistan", "Sri Lanka",
    # Central Sub-Saharan Africa
    "Angola", "Central African Republic", "Republic of the Congo", "Democratic Republic of the Congo", "Equatorial Guinea", "Gabon",
    "Congo, Dem. Rep.",
    "Congo, Rep.",
    # Eastern Sub-Saharan Africa
    "Burundi", "Comoros", "Djibouti", "Eritrea", "Ethiopia", "Kenya", "Madagascar", "Malawi", "Mozambique", "Rwanda",
    "Somalia", "South Sudan", "Tanzania", "Uganda", "Zambia", "Mauritius", "Seychelles",
    # Southern Sub-Saharan Africa
    "Botswana", "Lesotho", "Namibia", "South Africa", "Eswatini", "Zimbabwe", "Western Sahara",
    # Western Sub-Saharan Africa
    "Benin", "Burkina Faso", "Cameroon", "Cape Verde", "Chad", "Ivory Coast", "Gambia", "Ghana", "Guinea", "Guinea-Bissau",
    "Liberia", "Mali", "Mauritania", "Niger", "Nigeria", "Senegal", "Sierra Leone", "Togo", "São Tomé and Príncipe",
    "Cote d'Ivoire",  # Ivory Coast is also called "Cote d'Ivoire" in some contexts
    "Cote D'Ivoire",  # Ivory Coast is also called "Cote d'Ivoire" in some contexts
    "Gambia, The",
    "Sao Tome and Principe"
  ),
  
  ihme_region = c(
    # East Asia
    rep("Southeast Asia, East Asia, and Oceania", 1+6+12+21+2),
    
    # Central Europe
    rep("Central Europe, Eastern Europe, and Central Asia", 13+7+8+4),
    
    # Western Europe
    rep("High-income", 23+2),
    
    # Caribbean
    rep("Latin America and Caribbean", 18+3+10+3+2+2),
    
    # North Africa and Middle East
    rep("North Africa and Middle East", 22),
    
    # South Asia
    rep("South Asia", 6),
    
    # Central Sub-Saharan Africa
    rep("Sub-Saharan Africa", 6+17+9+19+4)
  ),
  lower_ihme_region = c(
    # East Asia 7+2
    rep("East Asia", 1+3), # Hong Kong, China, Taiwan, North Korea,
    rep("High-income Asia Pacific", 2), # Japan, South Korea,
    rep("East Asia", 1+2), # Mongolia, Korea Rep., South Korea

    # Southeast Asia
    rep("Southeast Asia", 6), # Cambodia through Phillipines
    rep("High-income Asia Pacific", 1), # Singapore
    rep("Southeast Asia", 3), # Thailand, Timor-Leste, Vietnam
    rep("High-income Asia Pacific", 1), # Brunei
    rep("South Asia", 1), # Maldives

    # Oceania
    rep("Oceania", 1), # American Samoa
    rep("Australasia", 1), # Australia
    rep("Oceania", 7), # Fiji
    rep("Oceania", 1), # New Zealand
    rep("Oceania", 11), # Fiji through French Polynesia

    # Central Asia
    rep("Central Asia", 8+1),
    
    # Central Europe
    rep("Central Europe", 13+3),

    # Eastern Europe
    rep("Eastern Europe", 7),
    
    # Western Europe
    rep("Western Europe", 23),

    # High-income North America
    rep("High-income North America", 2),
    
    # Caribbean
    rep("Caribbean", 18+2),

    # Andean Latin America
    rep("Andean Latin America", 3),

    # Central Latin America
    rep("Central Latin America", 10),

    # Southern Latin America
    rep("Southern Latin America", 3),

    # Tropical Latin America
    rep("Tropical Latin America", 2),
    
    # North Africa and Middle East
    rep("North Africa and Middle East", 22),
    
    # South Asia
    rep("South Asia", 6),
    
    # Central Sub-Saharan Africa
    rep("Central Sub-Saharan Africa", 6+2),

    # Eastern Sub-Saharan Africa
    rep("Eastern Sub-Saharan Africa", 17),

    # Southern Sub-Saharan Africa
    rep("Southern Sub-Saharan Africa", 7),

    # Western Sub-Saharan Africa
    rep("Western Sub-Saharan Africa", 19+4)
  ),
   lending_group = c(
    # East Asia
    "High income",       # Hong Kong SAR, China
    "Upper middle income", # China
    "High income",       # Taiwan, China
    "Low income",        # Korea, North (DPR Korea - generally considered low income, though data can be limited)
    "High income",       # Japan
    "High income",       # Korea, Rep. (South Korea)
    "Upper middle income", # Mongolia
    "High income",       # Korea, Rep. (South Korea) - Duplicate
    "High income",       # South Korea - Duplicate
    # Southeast Asia
    "Lower middle income", # Cambodia
    "Lower middle income", # Indonesia
    "Lower middle income", # Lao PDR (Laos)
    "Upper middle income", # Malaysia
    "Lower middle income", # Myanmar
    "Lower middle income", # Philippines
    "High income",       # Singapore
    "Upper middle income", # Thailand
    "Lower middle income", # Timor-Leste
    "Lower middle income", # Vietnam
    "High income",       # Brunei Darussalam
    "Upper middle income", # Maldives
    # Oceania
    "Upper middle income", # American Samoa
    "High income",       # Australia
    "Upper middle income", # Fiji
    "Lower middle income", # Micronesia, Fed. Sts. (Micronesia)
    "High income",       # Guam
    "Lower middle income", # Kiribati
    "Upper middle income", # Marshall Islands
    "High income",       # Nauru
    "High income",       # New Caledonia
    "High income",       # New Zealand
    "Upper middle income", # Niue (Often grouped with other small island developing states, often upper-middle)
    "High income",       # Northern Mariana Islands
    "High income",       # Palau
    "Lower middle income", # Papua New Guinea
    "Lower middle income", # Samoa
    "Lower middle income", # Solomon Islands
    "Upper middle income", # Tonga
    "Upper middle income", # Tuvalu
    "Lower middle income", # Vanuatu
    "High income",       # Cook Islands (Often grouped with high-income Pacific islands)
    "High income",       # French Polynesia
    # Central Asia
    "Upper middle income", # Kazakhstan
    "Lower middle income", # Kyrgyz Republic (Kyrgyzstan)
    "Lower middle income", # Tajikistan
    "Upper middle income", # Turkmenistan
    "Lower middle income", # Uzbekistan
    "Upper middle income", # Georgia
    "Upper middle income", # Armenia
    "Upper middle income", # Azerbaijan
    "Lower middle income", # Kyrgyz Republic - Duplicate
    # Central Europe
    "Upper middle income", # Albania
    "High income",       # Bulgaria
    "Upper middle income", # Bosnia and Herzegovina
    "High income",       # Croatia
    "High income",       # Czechia (Czech Republic)
    "High income",       # Hungary
    "Upper middle income", # North Macedonia (Macedonia)
    "Upper middle income", # Montenegro
    "High income",       # Poland
    "Upper middle income", # Romania
    "Upper middle income", # Serbia
    "High income",       # Slovak Republic (Slovakia)
    "High income",       # Slovenia
    "High income",       # Czechia - Duplicate
    "High income",       # Slovak Republic - Duplicate
    "Upper middle income", # North Macedonia - Duplicate (Macedonia)
    # Eastern Europe
    "Upper middle income", # Belarus
    "High income",       # Estonia
    "High income",       # Latvia
    "High income",       # Lithuania
    "Lower middle income", # Moldova
    "High income",       # Russian Federation (Russia)
    "Lower middle income", # Ukraine
    # Western Europe
    "High income",       # Andorra
    "High income",       # Austria
    "High income",       # Belgium
    "High income",       # Cyprus
    "High income",       # Denmark
    "High income",       # Finland
    "High income",       # France
    "High income",       # Germany
    "High income",       # Greece
    "High income",       # Iceland
    "High income",       # Ireland
    "High income",       # Italy
    "High income",       # Luxembourg
    "High income",       # Malta
    "High income",       # Monaco
    "High income",       # Netherlands
    "High income",       # Norway
    "High income",       # Portugal
    "High income",       # San Marino
    "High income",       # Spain
    "High income",       # Sweden
    "High income",       # Switzerland
    "High income",       # United Kingdom
    # High-income North America
    "High income",       # Canada
    "High income",       # United States
    # Caribbean
    "High income",       # Antigua and Barbuda
    "High income",       # Bahamas, The (Bahamas)
    "High income",       # Barbados
    "Upper middle income", # Cuba (Data limited, but generally placed here)
    "Upper middle income", # Dominican Republic
    "Upper middle income", # Grenada
    "Low income",        # Haiti
    "Upper middle income", # Jamaica
    "High income",       # Puerto Rico
    "Upper middle income", # St. Lucia
    "Upper middle income", # St. Vincent and the Grenadines
    "High income",       # Trinidad and Tobago
    "High income",       # British Virgin Islands
    "High income",       # US Virgin Islands
    "Upper middle income", # Dominica
    "High income",       # St. Kitts and Nevis
    "Upper middle income", # Guyana
    "Upper middle income", # Suriname
    "Upper middle income", # St. Lucia - Duplicate
    "Upper middle income", # St. Vincent and the Grenadines - Duplicate
    # Andean Latin America
    "Lower middle income", # Bolivia
    "Upper middle income", # Ecuador
    "Upper middle income", # Peru
    # Central Latin America
    "Upper middle income", # Colombia
    "Upper middle income", # Costa Rica
    "Lower middle income", # El Salvador
    "Lower middle income", # Guatemala
    "Lower middle income", # Honduras
    "Upper middle income", # Mexico
    "Lower middle income", # Nicaragua
    "High income",       # Panama
    "Upper middle income", # Venezuela, RB (Venezuela) - *Note: World Bank currently lists Venezuela as 'Unclassified' due to data unavailability, but historically upper-middle income.*
    "Upper middle income", # Belize
    # Southern Latin America
    "Upper middle income", # Argentina
    "High income",       # Chile
    "High income",       # Uruguay
    # Tropical Latin America
    "Upper middle income", # Brazil
    "Upper middle income", # Paraguay
    # North Africa and Middle East
    "Low income",        # Afghanistan
    "Lower middle income", # Algeria
    "High income",       # Bahrain
    "Lower middle income", # Egypt, Arab Rep. (Egypt)
    "Upper middle income", # Iran, Islamic Rep. (Iran)
    "Upper middle income", # Iraq
    "Lower middle income", # Jordan
    "High income",       # Kuwait
    "Lower middle income", # Lebanon
    "Upper middle income", # Libya
    "Lower middle income", # Morocco
    "High income",       # Oman
    "Lower middle income", # West Bank and Gaza (Palestine)
    "High income",       # Qatar
    "High income",       # Saudi Arabia
    "Low income",        # Sudan
    "Low income",        # Syrian Arab Republic (Syria)
    "Lower middle income", # Tunisia
    "Upper middle income", # Türkiye (Turkey)
    "High income",       # United Arab Emirates
    "Low income",        # Yemen
    "High income",       # Israel
    # South Asia
    "Lower middle income", # Bangladesh
    "Lower middle income", # Bhutan
    "Lower middle income", # India
    "Lower middle income", # Nepal
    "Lower middle income", # Pakistan
    "Lower middle income", # Sri Lanka
    # Central Sub-Saharan Africa
    "Lower middle income", # Angola
    "Low income",        # Central African Republic
    "Lower middle income", # Congo, Rep. (Republic of the Congo)
    "Low income",        # Congo, Dem. Rep. (Democratic Republic of the Congo)
    "Upper middle income", # Equatorial Guinea
    "Upper middle income", # Gabon
    "Low income",        # Congo, Dem. Rep. - Duplicate
    "Lower middle income", # Congo, Rep. - Duplicate
    # Eastern Sub-Saharan Africa
    "Low income",        # Burundi
    "Lower middle income", # Comoros
    "Lower middle income", # Djibouti
    "Low income",        # Eritrea
    "Low income",        # Ethiopia
    "Lower middle income", # Kenya
    "Low income",        # Madagascar
    "Low income",        # Malawi
    "Low income",        # Mozambique
    "Low income",        # Rwanda
    "Low income",        # Somalia
    "Low income",        # South Sudan
    "Low income",        # Tanzania
    "Low income",        # Uganda
    "Low income",        # Zambia
    "High income",       # Mauritius
    "High income",       # Seychelles
    # Southern Sub-Saharan Africa
    "Upper middle income", # Botswana
    "Lower middle income", # Lesotho
    "Upper middle income", # Namibia
    "Upper middle income", # South Africa
    "Lower middle income", # Eswatini
    "Low income",        # Zimbabwe
    "Unclassified",      # Western Sahara (Not classified by World Bank due to disputed status)
    # Western Sub-Saharan Africa
    "Lower middle income", # Benin
    "Low income",        # Burkina Faso
    "Lower middle income", # Cameroon
    "Lower middle income", # Cabo Verde (Cape Verde)
    "Low income",        # Chad
    "Lower middle income", # Côte d'Ivoire (Ivory Coast)
    "Low income",        # Gambia, The (Gambia)
    "Lower middle income", # Ghana
    "Low income",        # Guinea
    "Low income",        # Guinea-Bissau
    "Low income",        # Liberia
    "Low income",        # Mali
    "Lower middle income", # Mauritania
    "Low income",        # Niger
    "Lower middle income", # Nigeria
    "Lower middle income", # Senegal
    "Low income",        # Sierra Leone
    "Low income",        # Togo
    "Lower middle income", # Sao Tome and Principe
    "Lower middle income", # Cote d'Ivoire - Duplicate
    "Lower middle income", # Cote D'Ivoire - Duplicate
    "Low income",        # Gambia, The - Duplicate
    "Lower middle income"  # Sao Tome and Principe - Duplicate
)

)

# Display the first few rows of the mapping
head(iso3_ihme_mapping)

# You can save this mapping to a CSV file if needed
# write.csv(iso3_ihme_mapping, "iso3_ihme_mapping.csv", row.names = FALSE)

# Sample usage: Find the IHME region for a specific ISO3 code
get_ihme_region <- function(iso3_code) {
  region <- iso3_ihme_mapping$ihme_region[iso3_ihme_mapping$iso3 == iso3_code]
  if (length(region) == 0) {
    return(NA)
  } else {
    return(region)
  }
}

atc_mapping <- list(
  "J01A" = c("Tetracycline", "Tetracyclines", "Minocycline", "Tigecycline", "Glycylcyclines"),
  
  "J01C" = c("Beta lactam antibacterials, penicillins", "Penicillins", "BSP", "NSP", "Mecillinam",
             "Broad spectrum penicillins",  "Narrow spectrum penicillins",
             "Aminopenicillin", "Anti-pseudomonal penicillin/Beta-Lactamase inhibitors",
             "Beta Lactam/Beta-lactamase inhibitors", "Penicillin", "Methicillin", 
             "Amoxy/clav", "Ampicillin", "Oxacillin", "Ampicillin sulbactam", "Pip/taz"),
  
  "J01D" = c("Cephalosporins", "Cephalosporin (undefined and 1st/2nd gen)", "Cephalosporins (3rd/4th gen)",
             "Third-generation cephalosporins", "Fourth-generation cephalosporins",
             "Carbapenems", "Imipenem", "Meropenem", "Ceftriaxone", "Cefotaxime", 
             "Ceftazidime", "Cefixime", "Doripenem", "Ceftazidime avibactam", 
             "Cefepime", "Aztreonam", "Ceftolozane tazobactam", "Meropenem vaborbactam",
             "Cefpodoxime", "Ceftaroline", "Ceftibuten", "Ertapenem", "Monobactams"),
  
  "J01E" = c("Sulfonamide and Trimethoprim", "Sulfamethoxaole", "Trimethoprim",
             "Trimethoprim-Sulfamethoxazole", "Trimethoprim sulfa"),
  
  "J01F" = c("Macrolides, Lincosamides and Streptogramins", "Macrolide", "Macrolides", "Azithromycin", 
             "Clarithromycin", "Erythromycin", "Telithromycin", "Clindamycin", 
             "Quinupristin dalfopristin"),
  
  "J01G" = c("Aminoglycosides", "Gentamicin", "Amikacin", "Tobramycin", "Neomycin", 
             "Streptomycin"),
  
  "J01M" = c("Fluoroquinolones", "Quinolones", "Nalidixic Acid", "Ciprofloxacin", 
             "Levofloxacin", "Ofloxacin", "Moxifloxacin"),
  
  "J01X" = c("Nitrofuran derivatives", "Vancomycin", "Glycopeptides", "Polymyxins", 
             "Daptomycin", "Linezolid", "Colistin", "Teicoplanin", "Lipopeptides","Phosphonics"),
  
  "Other" = c("MDR", "J01XX Other antibacterials", "Resistance to one or more antibiotics", 
              "Rifampin", "Oxazolidinones")
)
get_atc_class <- function(atc_code) {
  for (class in names(atc_mapping)) {
    if (atc_code %in% atc_mapping[[class]]) {
      return(class)
    }
  }
  return(NA)
}
atc_names <- list(
  "J01A" = "Tetracyclines",
  "J01B" = "Glycopeptides and Lipopeptides",
  "J01C" = "Penicillins",
  "J01D" = "Other Beta-Lactams",
  "J01E" = "Sulfonamides and Trimethoprim",
  "J01F" = "Macrolides, Lincosamides and Streptogramins",
  "J01G" = "Aminoglycosides",
  "J01M" = "Quniolones"
)

bacteria_mapping <- data.frame(
  in_names = c(
    "Acinetobacter baumannii", # Maps to Acinetobacter spp.
    "Escherichia coli", # Maps to E. coli
    "Klebsiella pneumoniae", # Maps to K. pneumoniae
    "Streptococcus pneumoniae", # Maps to S. pneumoniae
    "Non-typhoidal Salmonella", # Maps to Salmonella spp.
    "Salmonella Typhi", # Maps to Salmonella spp.
    "Salmonella Paratyphi", # Maps to Salmonella spp.
    "Shigella spp.", # Maps to Shigella spp.
    "Staphylococcus aureus", # Maps to S. aureus
    "Enterococcus faecalis", # Maps to E. faecalis
    "Enterococcus faecium", # Maps to E. faecium
    "Pseudomonas aeruginosa", # Maps to P. aeruginosa
    "Enterobacter spp.",
    "Group A Streptococcus", # Maps to S. pyogenes
    "Group B Streptococcus", # Maps to S. agalactiae
    "Other enterococci", # Maps to Enterococcus spp.
    "Proteus spp.",
    "Citrobacter spp.",
    "Serratia spp.",
    "Haemophilus influenzae", # Maps to H. influenzae
    "Morganella spp.",
    "Streptococcus pyogenes", # Maps to S. pyogenes
    "Streptococcus agalactiae", # Maps to S. agalactiae
    "Shigella spp", # Duplicate without period
    "Enterobacter spp", # Duplicate without period
    "Proteus spp", # Duplicate without period
    "Citrobacter spp", # Duplicate without period
    "Serratia spp", # Duplicate without period
    "Morganella spp", # Duplicate without period
    "Enterococcus spp", # Maps to Enterococcus spp.
    "Salmonella spp", # Maps to Salmonella spp.
    # Additional mappings based on your frequency table
    "A. baumannii", # Alternative name for Acinetobacter baumannii
    "Acinetobacter", # Maps to Acinetobacter spp.
    "Enterobacteriaceae", # Family level, maps to itself
    "Enterococci", # Maps to Enterococcus spp.
    "Gram-negatives", # Broad category, maps to itself
    "H. influenza", # Alternative spelling, maps to H. influenzae
    "Klebsiella spp.", # Maps to K. pneumoniae (most common species)
    "MRSA", # Methicillin-resistant S. aureus, maps to S. aureus
    "Mycobacterium tuberculosis", # Maps to itself
    "N. meningitis", # Maps to N. meningitidis
    "P. mirabilis", # Species of Proteus, maps to Proteus spp.
    "Pseudomonas spp.", # Maps to P. aeruginosa (most common species)
    "S. epidermydis", # Maps to itself
    "Strep. Viridians", # Maps to Viridans group streptococci
    "Morganella morganii",
    "Enterococcus non-speciated",
    "Neisseria gonorrhoeae",
    "Enterococcus spp."
  ),
  joe_names = c(
    "Acinetobacter spp.",
    "E. coli",
    "K. pneumoniae",
    "S. pneumoniae",
    "Salmonella spp.",
    "Salmonella spp.",
    "Salmonella spp.",
    "Shigella spp.",
    "S. aureus",
    "E. faecalis",
    "E. faecium",
    "P. aeruginosa",
    "Enterobacter spp.",
    "S. pyogenes",
    "S. agalactiae",
    "Enterococcus spp.",
    "Proteus spp.",
    "Citrobacter spp.",
    "Serratia spp.",
    "H. influenzae",
    "Morganella spp.",
    "S. pyogenes",
    "S. agalactiae",
    "Shigella spp.",
    "Enterobacter spp.",
    "Proteus spp.",
    "Citrobacter spp.",
    "Serratia spp.",
    "Morganella spp.",
    "Enterococcus spp.",
    "Salmonella spp.",
    # Additional mappings
    "Acinetobacter spp.",
    "Acinetobacter spp.",
    "Enterobacteriaceae",
    "Enterococcus spp.",
    "Gram-negatives",
    "H. influenzae",
    "K. pneumoniae",
    "S. aureus",
    "M. tuberculosis",
    "N. meningitidis",
    "Proteus spp.",
    "P. aeruginosa",
    "S. epidermidis",
    "Streptococcus spp.",
    "Morganella spp.",
    "Enterococcus spp.",
    "N. gonorrhoeae",
    "Enterococcus spp."
  )
)
get_bacteria_name <- function(ihme_name) {
  if (ihme_name %in% bacteria_mapping$joe_names) {
    return(ihme_name)
  }
  match <- bacteria_mapping[bacteria_mapping$in_names == ihme_name, "joe_names"]
  if (length(match) > 0 && !is.na(match)) {
    return(match)
  } else {
    return(NA)
  }
}


