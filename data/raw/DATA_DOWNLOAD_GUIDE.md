# Data Download Guide

This document provides step-by-step instructions for downloading each dataset needed for the Portugal Earthquake project. All files should be placed in their corresponding subdirectory under `data/raw/`.

---

## 1. Slave Voyages — Trans-Atlantic Slave Trade Database

**What it contains:** Voyage-level data (36,000+ voyages) including ship nationality/flag, embarkation/disembarkation ports, number of enslaved people, dates, and 274 variables total. Key for measuring Portuguese vs. non-Portuguese commercial activity in Brazil over time.

**Where:** https://www.slavevoyages.org/voyage/downloads

**Steps:**
1. Go to https://www.slavevoyages.org/voyage/downloads
2. Download the full Trans-Atlantic database in CSV format
3. Save to: `data/raw/slave_voyages/`

**Alternative (filtered download):**
1. Go to https://www.slavevoyages.org/voyage/database/
2. Use filters to select Portuguese-flagged voyages, or voyages to/from Brazilian ports
3. Click Download > Excel > Filtered results with visible columns
4. Save to: `data/raw/slave_voyages/`

**Key variables for this project:**
- `NATINIMP` — Flag/nationality of vessel (filter for Portuguese)
- `MJBYPTIMP` — Principal place of slave purchase
- `MJSLPTIMP` — Principal place of slave landing
- `YEARAM` — Year of arrival at port of disembarkation
- `SLAXIMP` — Total slaves embarked (imputed)
- `SLAMIMP` — Total slaves disembarked (imputed)

---

## 2. TePaske/Klein — Royal Treasuries of the Spanish Empire

**What it contains:** Annual fiscal data (revenue by tax type, expenditures) for individual Spanish colonial treasury offices (cajas reales), ~1570s-1820s. This is the control group fiscal data.

**Where:** https://realhacienda.colmex.mx/

**Steps:**
1. Go to https://realhacienda.colmex.mx/
2. Navigate to the database/data section (look for "Base de datos" or the search tool)
3. The data is in Excel format — download files for the regions of interest:
   - New Spain (Mexico)
   - Peru
   - New Granada (Colombia)
   - Rio de la Plata (Argentina)
   - Chile
4. Save all files to: `data/raw/tepaske_klein/`

**Alternative source:** The original books are on Internet Archive:
https://archive.org/details/royaltreasurieso0002tepa

**Key variables:**
- Treasury office (caja) identifier
- Year
- Revenue category (tributo, alcabala, almojarifazgo, quinto real, diezmos, tobacco monopoly, etc.)
- Revenue amount
- Expenditure category and amount

---

## 3. Nuno Palma — Portuguese GDP and Economic Data

**What it contains:** Portuguese GDP series 1527-1850, plus prices, wages, and rents data from 1300-1910. Key for the first stage (measuring Portuguese state capacity decline).

### 3a. GDP Data
**Where:** https://sites.google.com/site/npgpalma/publications

**Steps:**
1. Go to Nuno Palma's publications page (link above)
2. Find the entry for "From Convergence to Divergence: Portuguese Economic Growth, 1527-1850" (JEH, 2019)
3. Download the replication files / data appendix (should include Excel with GDP series)
4. Save to: `data/raw/palma_portuguese_gdp/`

### 3b. Prices, Wages, Rents
**Where:** http://pwr-portugal.ics.ul.pt/

**Steps:**
1. Go to the PWR Portugal website (link above)
2. Navigate to "The PWR Data Files" section
3. Download the Excel files for prices, wages, and rents
4. Save to: `data/raw/palma_portuguese_gdp/`

**Key variables:**
- Real GDP per capita (annual, 1527-1850)
- Nominal GDP
- CPI / price indices
- Commodity prices (134 commodities across Lisbon, Porto, Coimbra, Evora)
- Wages (26 occupational categories)

---

## 4. Karaman & Pamuk — European State Finances

**What it contains:** Revenue per capita and revenue/GDP ratios for European states including Portugal and Spain, 1500-1914. Key for comparing Portuguese fiscal capacity to other European powers.

**Where:** Published as supplementary data to:
Karaman, K. Kivanc, and Sevket Pamuk. "Ottoman State Finances in European Perspective, 1500-1914." *Journal of Economic History* 70, no. 3 (2010): 593-629.

**Steps:**
1. Access the paper at Cambridge Core or the PDF at:
   http://piketty.pse.ens.fr/files/KaramanPamuk2010.pdf
2. The data appendix tables contain revenue/capita figures for:
   - Portugal, Spain, England, France, Netherlands, Ottoman Empire, and others
3. If a standalone dataset exists, check:
   https://ata.bogazici.edu.tr (Bogazici University — Pamuk's institution)
4. May need to hand-digitize tables from the paper/appendix
5. Save to: `data/raw/karaman_pamuk/`

**Key variables:**
- Country
- Year (benchmark years, not annual for all states)
- Revenue per capita (in grams of silver)
- Revenue/GDP ratio

---

## 5. Maddison Project Database 2023

**Status: ALREADY DOWNLOADED** -> `data/raw/maddison/mpd2023.xlsx`

**What it contains:** GDP per capita estimates for 169 countries, from antiquity to 2022. Includes Portugal, Spain, Brazil (from 1800).

**Source:** https://www.rug.nl/ggdc/historicaldevelopment/maddison/releases/maddison-project-database-2023

---

## 6. Carrara — Brazilian Colonial Fiscal Data

**What it contains:** Revenue and expenditure data for Brazilian captaincies during the colonial period (17th-18th century). This is the treatment group fiscal data — the Portuguese-side equivalent of TePaske/Klein.

**Where:** This data is published in book form and is NOT available as a digital download.

**Primary source:**
Carrara, Angelo Alves. *Receitas e Despesas da Real Fazenda no Brasil, Seculo XVIII.* Juiz de Fora: Editora UFJF, 2009.

**Steps:**
1. Obtain the book via library loan or purchase
2. Hand-digitize the fiscal tables into a spreadsheet
3. Alternatively, contact Angelo Carrara at UFJF (carrara@pq.cnpq.br or check his faculty page at ufjf.br) to request the underlying data
4. Save digitized data to: `data/raw/carrara_brazil_fiscal/`

**Related published article with summary data:**
Carrara, Angelo Alves. "Costos y beneficios de una colonia: introduccion a la fiscalidad colonial del Estado de Brasil, 1607-1808." *Investigaciones de Historia Economica* 6 (2010).
https://www.sciencedirect.com/science/article/abs/pii/S1698698910700333

**Key variables:**
- Captaincy (administrative unit)
- Year
- Revenue by source (dizimos/tithes, entradas/transit duties, quinto do ouro/gold royalty, monopoly revenues)
- Expenditure by category
- Remittances to Lisbon

---

## 7. Additional Sources (Lower Priority)

### Naval Tonnage — Glete (1993)
- Source: Glete, Jan. *Navies and Nations: Warships, Navies and State Building in Europe and America, 1500-1860.* Stockholm: Almqvist & Wiksell, 1993.
- Contains: Decadal naval tonnage by European power (Portugal, Spain, Britain, France, Netherlands)
- Access: Library; hand-digitize appendix tables
- Save to: `data/raw/` (create subdirectory as needed)

### Gold Shipments — Noya Pinto (1979)
- Source: Noya Pinto, Virgilio. *O Ouro Brasileiro e o Comercio Anglo-Portugues.* Sao Paulo: Companhia Editora Nacional, 1979.
- Contains: Annual Brazilian gold shipments to Portugal (18th century)
- Access: Library; hand-digitize
- Save to: `data/raw/` (create subdirectory as needed)

### Fisher (1971) — Anglo-Portuguese Trade
- Source: Fisher, H.E.S. *The Portugal Trade: A Study of Anglo-Portuguese Commerce, 1700-1770.* London: Methuen, 1971.
- Contains: Annual trade volumes between England and Portugal
- Access: Library; hand-digitize tables

### Dincecco (2009) — Fiscal Centralization
- Source: Dincecco, Mark. "Fiscal Centralization, Limited Government, and Public Revenues in Europe, 1650-1913." *JEH* 69, no. 1 (2009).
- Contains: Fiscal centralization measures and revenue/capita for European states
- Access: Journal supplementary materials

### COLDAT — Colonial Dates Dataset
- Contains: Colonizer identity, colonial duration, independence dates
- Access: Search Harvard Dataverse for "COLDAT"

### ICOW — Colonial History Dataset
- Contains: Colonial claims, territorial changes
- Access: https://sites.psu.edu/icow/
