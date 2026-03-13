## ─────────────────────────────────────────────────────────────────────────────
## Visualise many-to-one & many-to-many mappings from a curation map CSV
## as polished bipartite diagrams.
##
## Usage:  Rscript plot_mapping_types.R <path/to/map.csv>
## ─────────────────────────────────────────────────────────────────────────────

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ── Locate data (CLI arg or default) ─────────────────────────────────────────
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1 && file.exists(cli_args[1])) {
    csv_path <- cli_args[1]
} else {
    script_dir <- tryCatch(
        dirname(rstudioapi::getSourceEditorContext()$path),
        error = function(e) {
            args <- commandArgs(trailingOnly = FALSE)
            file_arg <- grep("^--file=", args, value = TRUE)
            if (length(file_arg)) dirname(sub("^--file=", "", file_arg))
            else getwd()
        }
    )
    csv_path <- file.path(script_dir, "cBioPortal_treatment_type_map.csv")
    if (!file.exists(csv_path))
        csv_path <- "cBioPortalData/maps/cBioPortal_treatment_type_map.csv"
}

message("Reading → ", csv_path)
dat <- read.csv(csv_path, stringsAsFactors = FALSE)

# ── Derive a human-readable label from the filename ──────────────────────────
csv_basename <- tools::file_path_sans_ext(basename(csv_path))
pretty_label <- csv_basename |>
    sub("^cBioPortal_", "", x = _) |>
    sub("_map$", "", x = _) |>
    gsub("_", " ", x = _) |>
    tools::toTitleCase()

# ── Auto-detect the curated ontology column ──────────────────────────────────
col_names <- names(dat)
if ("curated_ontology" %in% col_names) {
    curated_col <- "curated_ontology"
} else {
    # R converts ":" to "." in column names; match both patterns
    candidates <- grep("[.:]curated_ontology$", col_names, value = TRUE)
    preferred  <- grep("^sample_type[.:]curated_ontology$", candidates,
                       value = TRUE)
    if (length(preferred) > 0) {
        curated_col <- preferred[1]
    } else if (length(candidates) > 0) {
        curated_col <- candidates[1]
    } else {
        stop("No curated_ontology column found in: ", csv_path)
    }
}
message("Using curated column: ", curated_col)

# Normalise into a two-column frame (original_value, curated_ontology)
dat_raw <- dat %>%
    select(original_value, curated_ontology = all_of(curated_col)) %>%
    filter(!is.na(curated_ontology), curated_ontology != "")

# ── Auto-detect the multi-value delimiter in curated_ontology ────────────────
vals <- dat_raw$curated_ontology
if (any(grepl("<;>", vals, fixed = TRUE))) {
    delim <- "<;>"
} else if (any(grepl(";", vals, fixed = TRUE))) {
    delim <- ";"
} else {
    delim <- NULL
}
if (!is.null(delim)) message("Detected delimiter: '", delim, "'")

# ── Expand delimited curated_ontology + deduplicate rows ─────────────────────
if (!is.null(delim)) {
    # Escape the delimiter for use as regex in separate_rows
    delim_regex <- gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", delim)
    dat_long <- dat_raw %>%
        mutate(row_id = row_number()) %>%
        separate_rows(curated_ontology, sep = delim_regex) %>%
        mutate(curated_ontology = trimws(curated_ontology)) %>%
        filter(curated_ontology != "") %>%
        select(original_value, curated_ontology) %>%
        distinct()
} else {
    dat_long <- dat_raw %>%
        mutate(curated_ontology = trimws(curated_ontology)) %>%
        select(original_value, curated_ontology) %>%
        distinct()
}

# Forward cardinality
fwd <- dat_long %>%
    group_by(original_value) %>%
    summarise(n_curated = n_distinct(curated_ontology), .groups = "drop")

# ── Many-to-one edges ───────────────────────────────────────────────────────
pure_one <- fwd %>% filter(n_curated == 1) %>% pull(original_value)

m2o_edges <- dat_long %>%
    filter(original_value %in% pure_one) %>%
    select(original_value, curated_ontology) %>%
    distinct()

m2o_groups <- m2o_edges %>%
    group_by(curated_ontology) %>%
    filter(n() > 1) %>%
    mutate(group_size = n()) %>%
    ungroup()

top3 <- m2o_groups %>%
    distinct(curated_ontology, group_size) %>%
    slice_max(group_size, n = 3) %>%
    pull(curated_ontology)

m2o_sel <- m2o_groups %>% filter(curated_ontology %in% top3)

# ── Many-to-many edges ──────────────────────────────────────────────────────
m2m_originals <- fwd %>% filter(n_curated > 1) %>% pull(original_value)

m2m_edges_all <- dat_long %>%
    filter(original_value %in% m2m_originals) %>%
    select(original_value, curated_ontology) %>%
    distinct()

m2m_originals_sel <- m2m_edges_all %>%
    group_by(original_value) %>%
    summarise(n = n(), .groups = "drop") %>%
    slice_max(n, n = 6) %>%
    pull(original_value)

m2m_sel <- m2m_edges_all %>% filter(original_value %in% m2m_originals_sel)

# ── Colour palettes ─────────────────────────────────────────────────────────
pal_left  <- "#2b7bba"
pal_right <- "#e85d04"

edge_pal <- c("#6a4c93", "#1982c4", "#8ac926", "#ff595e", "#ff924c", "#56cfe1",
              "#457b9d", "#e63946", "#2a9d8f", "#f4a261")

# ── Bipartite plotting function ──────────────────────────────────────────────
bipartite_plot <- function(edges, title,
                           colour_by = c("right", "left"),
                           left_title  = "original_value",
                           right_title = "curated_ontology") {

    colour_by <- match.arg(colour_by)

    left_labs  <- unique(edges$original_value)
    right_labs <- unique(edges$curated_ontology)
    n_left     <- length(left_labs)
    n_right    <- length(right_labs)

    # Node positions
    left_df <- data.frame(
        label = left_labs, x = 0,
        y = seq(1, 0, length.out = n_left),
        stringsAsFactors = FALSE
    )
    right_df <- data.frame(
        label = right_labs, x = 1,
        y = seq(1, 0, length.out = n_right),
        stringsAsFactors = FALSE
    )

    # Segment coordinates
    seg_df <- edges %>%
        left_join(left_df  %>% rename(x1 = x, y1 = y),
                  by = c("original_value"  = "label")) %>%
        left_join(right_df %>% rename(x2 = x, y2 = y),
                  by = c("curated_ontology" = "label"))

    # Colour grouping column
    if (colour_by == "right") {
        seg_df$colour_group <- seg_df$curated_ontology
    } else {
        seg_df$colour_group <- seg_df$original_value
    }
    colour_levels <- unique(seg_df$colour_group)
    colour_map <- setNames(
        rep_len(edge_pal, length(colour_levels)),
        colour_levels
    )

    ggplot() +
        ## ── curved edges ────────────────────────────────────────────────
        geom_curve(
            data = seg_df,
            aes(x = x1, y = y1, xend = x2, yend = y2,
                colour = colour_group),
            curvature = 0.25, linewidth = 0.8, alpha = 0.6,
            show.legend = FALSE
        ) +
        scale_colour_manual(values = colour_map) +

        ## ── left labels (blue rounded boxes) ────────────────────────────
        geom_label(
            data = left_df,
            aes(x = x - 0.02, y = y, label = label),
            hjust        = 1,
            size         = 2.4,
            fontface     = "bold",
            colour       = "white",
            fill         = pal_left,
            linewidth    = 0,
            label.r      = unit(0.25, "lines"),
            label.padding = unit(0.2, "lines")
        ) +

        ## ── right labels (orange rounded boxes) ─────────────────────────
        geom_label(
            data = right_df,
            aes(x = x + 0.02, y = y, label = label),
            hjust        = 0,
            size         = 2.4,
            fontface     = "bold",
            colour       = "white",
            fill         = pal_right,
            linewidth    = 0,
            label.r      = unit(0.25, "lines"),
            label.padding = unit(0.2, "lines")
        ) +

        ## ── column headers ──────────────────────────────────────────────
        annotate("text", x = 0, y = 1.10, label = left_title,
                 fontface = "bold", size = 3.5, colour = "#444444") +
        annotate("text", x = 1, y = 1.10, label = right_title,
                 fontface = "bold", size = 3.5, colour = "#444444") +

        ## ── title & theme ───────────────────────────────────────────────
        labs(title = title) +
        coord_cartesian(xlim = c(-0.75, 1.75), ylim = c(-0.06, 1.20),
                        clip = "off") +
        theme_void(base_family = "sans") +
        theme(
            plot.title      = element_text(face = "bold", hjust = 0.5,
                                           size = 15, margin = margin(b = 6)),
            plot.background = element_rect(fill = "#fafafa", colour = NA),
            plot.margin     = margin(12, 12, 12, 12)
        )
}

# ── Build plots ──────────────────────────────────────────────────────────────
plots <- list()
has_m2o <- nrow(m2o_sel) > 0
has_m2m <- nrow(m2m_sel) > 0

if (has_m2o) plots$p1 <- bipartite_plot(m2o_sel, "Many-to-One Mapping",  colour_by = "right")
if (has_m2m) plots$p2 <- bipartite_plot(m2m_sel, "Many-to-Many Mapping", colour_by = "left")

if (length(plots) == 0) {
    message("No many-to-one or many-to-many mappings found in this file.")
    quit(save = "no")
}

ann_theme <- theme(
    plot.title    = element_text(face = "bold", size = 18,
                                 hjust = 0.5, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 11, hjust = 0.5,
                                  colour = "#666666",
                                  margin = margin(b = 10)),
    plot.background = element_rect(fill = "white", colour = NA)
)

if (length(plots) == 2) {
    combined <- plots$p1 + plots$p2 +
        plot_annotation(
            title    = paste0(pretty_label, " Curation Mappings"),
            subtitle = "cBioPortal original values → curated ontology terms",
            theme    = ann_theme
        )
    fig_w <- 20; fig_h <- 10
} else {
    combined <- plots[[1]] +
        plot_annotation(
            title    = paste0(pretty_label, " Curation Mappings"),
            subtitle = "cBioPortal original values → curated ontology terms",
            theme    = ann_theme
        )
    fig_w <- 12; fig_h <- 10
}

out_name <- paste0(csv_basename, "_plot.png")
out_path <- file.path(dirname(csv_path), out_name)
ggsave(out_path, combined, width = fig_w, height = fig_h, dpi = 200, bg = "white")
message("Saved → ", out_path)
