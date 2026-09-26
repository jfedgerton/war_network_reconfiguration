# ============================================================
# Stylized Ottoman-Byzantine territorial map, 1243-1400
# 2x2 panel figure for Gallipoli earthquake case
# REVISED v3: historical-accuracy fixes
#   - Köse Dağ now a battle marker, not a city
#   - Empire of Trebizond added as a polygon in panel a
#   - Battle of Bapheus (1301) marker added to panel b
#   - Battle of Çirmen (1371) marker added to panel d
# ============================================================
library(sf)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthdata)

# ------------------------------------------------------------
# Basemap: single land mass (no internal modern borders)
# ------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")
map_window <- st_as_sfc(
  st_bbox(
    c(xmin = 18, xmax = 44, ymin = 34, ymax = 48),
    crs = st_crs(4326)
  )
)
world_crop <- st_crop(world, map_window)
land_mask  <- st_make_valid(st_union(world_crop))   # one polygon, no inner borders

# ------------------------------------------------------------
# Period labels
period_1243 <- "1243: Seljuk defeat and Nicaean retreat"
period_1300 <- "c. 1300: Ottoman frontier polity emerges"
period_1352 <- "1354: Pre-earthquake battle lines"
period_1400 <- "1400: Pre-Ankara Ottoman peak"

# ------------------------------------------------------------
# Approximate territorial polygons
# Consolidated polity labels (5 categories instead of 13)
# ------------------------------------------------------------
territories <- st_sf(
  period = c(
    period_1243, period_1243, period_1243, period_1243, period_1243,
    period_1300, period_1300, period_1300, period_1300, period_1300,
    period_1352, period_1352, period_1352, period_1352,
    period_1400, period_1400, period_1400, period_1400, period_1400
  ),
  polity = c(
    "Byzantine state", "Seljuk Rum", "Byzantine state", "Balkan states", "Balkan states",
    "Byzantine state", "Ottoman state", "Turkmen / beyliks", "Balkan states", "Balkan states",
    "Byzantine state", "Ottoman state", "Balkan states", "Balkan states",
    "Byzantine state", "Byzantine state",
    "Ottoman state", "Ottoman state", "Balkan states"
  ),
  geometry = st_sfc(
    
    # --- 1243 Empire of Nicaea (NW Anatolia) ---
    st_polygon(list(matrix(c(
      26.5, 37.2,
      27.0, 39.5,
      28.5, 41.0,
      30.5, 41.0,
      31.0, 39.8,
      30.0, 37.8,
      28.0, 37.2,
      26.5, 37.2
    ), ncol = 2, byrow = TRUE))),
    
    # --- 1243 Seljuk Rum ---
    st_polygon(list(matrix(c(
      29.5, 36.5,
      38.5, 36.4,
      41.3, 39.0,
      39.0, 41.2,
      32.0, 41.0,
      29.5, 39.2,
      29.5, 36.5
    ), ncol = 2, byrow = TRUE))),
    
    # --- 1243 Empire of Trebizond (SE Black Sea coast) ---
    # Komnenian successor state, founded 1204; coloured as Byzantine
    st_polygon(list(matrix(c(
      38.5, 40.7,
      41.0, 40.7,
      41.4, 41.3,
      38.8, 41.3,
      38.5, 40.7
    ), ncol = 2, byrow = TRUE))),

    # --- 1243 Second Bulgarian Empire (dominant Balkan power under Ivan Asen II) ---
    st_polygon(list(matrix(c(
      22.5, 41.7,
      28.0, 41.9,
      28.0, 43.9,
      24.0, 44.2,
      22.0, 43.2,
      22.5, 41.7
    ), ncol = 2, byrow = TRUE))),

    # --- 1243 Western Balkan states (Serbia, Epirus / Thessalonica) ---
    st_polygon(list(matrix(c(
      19.0, 39.3,
      21.5, 39.5,
      22.0, 41.5,
      21.5, 43.3,
      19.5, 43.5,
      18.8, 41.5,
      19.0, 39.3
    ), ncol = 2, byrow = TRUE))),

    # --- 1300 Byzantine Empire (incl. Bithynia) ---
    st_polygon(list(matrix(c(
      22.0, 40.0,
      29.5, 40.0,
      30.5, 40.5,
      30.0, 41.2,
      28.5, 41.4,
      27.5, 41.7,
      23.0, 41.7,
      22.0, 41.0,
      22.0, 40.0
    ), ncol = 2, byrow = TRUE))),
    
    # --- 1300 Ottoman beylik (tight around Sogut/Bilecik) ---
    st_polygon(list(matrix(c(
      29.8, 39.8,
      30.6, 39.8,
      30.7, 40.4,
      30.0, 40.4,
      29.8, 39.8
    ), ncol = 2, byrow = TRUE))),
    
    # --- 1300 Other Anatolian beyliks ---
    st_polygon(list(matrix(c(
      26.0, 36.0,
      32.2, 36.2,
      32.0, 39.4,
      30.1, 40.0,
      27.0, 39.0,
      26.0, 36.0
    ), ncol = 2, byrow = TRUE))),

    # --- 1300 Second Bulgarian Empire (NE Balkans) ---
    st_polygon(list(matrix(c(
      23.0, 41.9,
      28.0, 42.1,
      28.0, 43.9,
      24.5, 44.2,
      23.0, 43.0,
      23.0, 41.9
    ), ncol = 2, byrow = TRUE))),

    # --- 1300 Kingdom of Serbia (NW Balkans, under Milutin) ---
    st_polygon(list(matrix(c(
      18.8, 41.3,
      22.5, 41.6,
      23.0, 43.3,
      21.0, 44.0,
      18.8, 43.3,
      18.8, 41.3
    ), ncol = 2, byrow = TRUE))),
    
    # --- 1352 Byzantine Empire ---
    st_polygon(list(matrix(c(
      22.8, 40.2,
      27.3, 40.1,
      29.1, 41.1,
      28.0, 41.9,
      24.0, 41.7,
      22.8, 40.2
    ), ncol = 2, byrow = TRUE))),
    
    # --- 1352 Ottoman polity ---
    st_polygon(list(matrix(c(
      26.5, 39.4,
      31.5, 39.4,
      31.5, 40.8,
      29.5, 41.2,
      26.5, 40.5,
      26.5, 39.4
    ), ncol = 2, byrow = TRUE))),
    
    # --- 1352 Serbian Empire (Dusan) ---
    st_polygon(list(matrix(c(
      18.8, 39.8,
      22.5, 39.8,
      23.5, 42.5,
      21.0, 43.5,
      18.8, 42.5,
      18.8, 39.8
    ), ncol = 2, byrow = TRUE))),
    
    # --- 1352 Bulgarian Empire (Ivan Alexander) ---
    # Southern frontier raised to sit clearly north of Byzantine
    # Adrianople (lat 41.68); Adrianople stays Byzantine until the
    # Ottoman capture, c.1361-69.
    st_polygon(list(matrix(c(
      22.5, 42.2,
      28.0, 42.5,
      28.0, 44.5,
      23.0, 44.2,
      22.5, 42.2
    ), ncol = 2, byrow = TRUE))),
    
    # --- 1400 Byzantine remnants (Constantinople area) ---
    st_polygon(list(matrix(c(
      28.55, 40.75,
      29.25, 40.75,
      29.25, 41.20,
      28.55, 41.20,
      28.55, 40.75
    ), ncol = 2, byrow = TRUE))),

    # --- 1400 Despotate of Morea (Byzantine — Peloponnese) ---
    st_polygon(list(matrix(c(
      21.4, 36.7,
      23.2, 36.5,
      23.4, 37.5,
      22.8, 38.3,
      21.6, 38.2,
      21.4, 36.7
    ), ncol = 2, byrow = TRUE))),

    # --- 1400 Ottoman Empire: Balkans / Thrace ---
    st_polygon(list(matrix(c(
      20.0, 40.1,
      29.5, 40.1,
      29.7, 44.7,
      25.5, 45.2,
      22.0, 44.0,
      20.0, 42.0,
      20.0, 40.1
    ), ncol = 2, byrow = TRUE))),
    
    # --- 1400 Ottoman Empire: Anatolia ---
    st_polygon(list(matrix(c(
      26.2, 36.0,
      38.0, 36.3,
      38.5, 40.8,
      29.0, 41.4,
      26.2, 39.5,
      26.2, 36.0
    ), ncol = 2, byrow = TRUE))),
    
    # --- 1400 Balkan vassal / influence zone ---
    st_polygon(list(matrix(c(
      19.0, 42.0,
      26.0, 42.3,
      26.5, 45.5,
      22.0, 46.0,
      19.0, 44.5,
      19.0, 42.0
    ), ncol = 2, byrow = TRUE))),
    
    crs = 4326
  )
)

# Clip territories to land mass (no polygons over water)
territories <- st_make_valid(territories)
territories_land <- suppressWarnings(
  st_intersection(territories, land_mask)
)

# Force panel ordering
territories_land$period <- factor(
  territories_land$period,
  levels = c(period_1243, period_1300, period_1352, period_1400)
)

# ------------------------------------------------------------
# Hubs (cities)
# ------------------------------------------------------------
hubs <- tibble(
  period = c(
    period_1243, period_1243, period_1243, period_1243,
    period_1300, period_1300, period_1300, period_1300, period_1300,
    period_1352, period_1352, period_1352, period_1352, period_1352,
    period_1400, period_1400, period_1400, period_1400, period_1400,
    period_1400, period_1400
  ),
  hub = c(
    "Constantinople (Latin)", "Nikaia", "Konya", "Trebizond",
    "Constantinople", "Bursa", "Nikaia", "Nikomedia", "Söğüt",
    "Constantinople", "Gallipoli", "Tzympe", "Bursa", "Adrianople",
    "Constantinople", "Gallipoli", "Edirne", "Sofia", "Nicopolis",
    "Bursa", "Mystras"
  ),
  lon = c(
    28.98, 29.72, 32.49, 39.72,
    28.98, 29.06, 29.72, 29.92, 30.18,
    28.98, 26.67, 26.78, 29.06, 26.56,
    28.98, 26.67, 26.56, 23.32, 24.90,
    29.06, 22.37
  ),
  lat = c(
    41.01, 40.43, 37.87, 41.00,
    41.01, 40.19, 40.43, 40.77, 40.02,
    41.01, 40.41, 40.52, 40.19, 41.68,
    41.01, 40.41, 41.68, 42.70, 43.70,
    40.19, 37.07
  )
)
hubs$period <- factor(
  hubs$period,
  levels = c(period_1243, period_1300, period_1352, period_1400)
)

# ------------------------------------------------------------
# Battle markers (separate from cities; rendered with an X glyph)
#   - Köse Dağ 1243: Mongol defeat of Seljuks
#   - Bapheus 1301: first major Ottoman victory over Byzantines
#   - Çirmen (Maritsa) 1371: Serbian collapse, Balkan vassalage
# ------------------------------------------------------------
battles <- tibble(
  period = c(period_1243, period_1300, period_1400),
  battle = c("Battle of Köse Dağ (1243)",
             "Battle of Bapheus (1301)",
             "Battle of Çirmen (1371)"),
  # Markers sit at the true (approximate) battle sites; label overlap is
  # handled by geom_text_repel, not by displacing the points:
  #   Köse Dağ  ~ near Suşehri (E Anatolia)
  #   Bapheus   ~ on the plain near Nicomedia (İzmit)
  #   Çirmen    ~ Chernomen/Ormenio on the Maritsa, SW of Edirne
  lon = c(38.30, 30.10, 26.40),
  lat = c(40.05, 40.70, 41.45)
)
battles$period <- factor(
  battles$period,
  levels = c(period_1243, period_1300, period_1352, period_1400)
)

# ------------------------------------------------------------
# Earthquake marker (1352 panel only)
# ------------------------------------------------------------
earthquake_marker <- tibble(
  period = period_1352,
  lon = 26.67,
  lat = 40.41,
  label = "1354 earthquake"
)
earthquake_marker$period <- factor(
  earthquake_marker$period,
  levels = c(period_1243, period_1300, period_1352, period_1400)
)

# ------------------------------------------------------------
# Ottoman expansion arrows
# ------------------------------------------------------------
expansion_arrows <- tibble(
  period = c(
    period_1352,
    period_1400,
    period_1400,
    period_1400
  ),
  x    = c(29.05, 26.67, 26.56, 26.56),
  y    = c(40.20, 40.41, 41.68, 41.68),
  xend = c(26.78, 26.56, 23.32, 24.90),
  yend = c(40.36, 41.68, 42.70, 43.70),
  label = c(
    "Bridgehead near Gallipoli",
    "Move into Thrace",
    "Expansion toward Bulgaria",
    "Expansion toward Danube"
  )
)
expansion_arrows$period <- factor(
  expansion_arrows$period,
  levels = c(period_1243, period_1300, period_1352, period_1400)
)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot() +
  # land mass — coastline only, no modern country borders
  geom_sf(
    data = land_mask,
    fill = "grey94",
    color = "grey55",
    linewidth = 0.25
  ) +
  # territorial polygons (already clipped to land)
  geom_sf(
    data = territories_land,
    aes(fill = polity),
    color = "white",
    linewidth = 0.3,
    alpha = 0.80
  ) +
  geom_curve(
    data = expansion_arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    curvature = 0.25,
    linewidth = 0.45,
    arrow = arrow(length = unit(0.12, "inches")),
    inherit.aes = FALSE
  ) +
  geom_point(
    data = hubs,
    aes(x = lon, y = lat),
    size = 1.3,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = hubs,
    aes(x = lon, y = lat, label = hub),
    size = 4.5,
    min.segment.length = 0,
    segment.linewidth = 0.2,
    box.padding = 1.5,
    point.padding = 0.5,
    force = 9,
    force_pull = 0.3,
    max.overlaps = 40,
    max.iter = 8000,
    seed = 123,
    inherit.aes = FALSE
  ) +
  # earthquake marker on the 1352 panel
  geom_point(
    data = earthquake_marker,
    aes(x = lon, y = lat),
    shape = 4, size = 4, color = "black", stroke = 0.8,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = earthquake_marker,
    aes(x = lon, y = lat, label = label),
    size = 5.5, color = "red", fontface = "italic",
    nudge_y = -2.1, nudge_x = -2.25,
    segment.color = "red", segment.linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  # battle markers: X glyph + italic label
  geom_point(
    data = battles,
    aes(x = lon, y = lat),
    shape = 4, size = 4, color = "grey15", stroke = 0.9,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = battles,
    aes(x = lon, y = lat, label = battle),
    size = 4.25, color = "grey15", fontface = "italic",
    min.segment.length = 0,
    segment.linewidth = 0.2,
    box.padding = 1.3,
    point.padding = 0.5,
    force = 8,
    force_pull = 0.3,
    max.overlaps = 40,
    max.iter = 8000,
    seed = 1235,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ period, ncol = 2) +
  coord_sf(
    xlim = c(18, 44),
    ylim = c(34, 48),
    expand = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Byzantine state"    = "#5B8FD6",
      "Ottoman state"      = "#6F9F45",
      "Seljuk Rum"         = "#B8A06A",
      "Turkmen / beyliks"  = "orange",
      "Balkan states"      = "#A283C6"
    ),
    breaks = c(
      "Byzantine state",
      "Ottoman state",
      "Balkan states",
      "Seljuk Rum",
      "Turkmen / beyliks"
    )
  ) +
  labs(
    #title = "Ottoman Expansion and Byzantine Retrenchment around the Gallipoli Earthquake",
    #subtitle = "Stylized territorial schematic, not a precise historical boundary reconstruction",
    fill = NULL,
    x = NULL,
    y = NULL,
    caption = paste(
      "Stylized territorial schematic, not a precise historical boundary reconstruction",
      sep = " "
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    strip.text = element_text(face = "bold", size = 18),
    strip.background = element_rect(fill = "grey97", color = "grey80"),
    plot.title = element_text(face = "bold", size = 24),
    plot.caption = element_text(size = 15, hjust = 0, color = "grey30"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = 12)
  )

print(p)

ggsave(
  filename = "results/main/ottoman_byzantine_four_panel_map.png",
  plot = p,
  width = 14, height = 11, dpi = 800
)
ggsave(
  filename = "results/main/ottoman_byzantine_four_panel_map.pdf",
  plot = p,
  width = 14, height = 12
)

