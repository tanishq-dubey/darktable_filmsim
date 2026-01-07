/*
 * Film Simulation Plugin for darktable
 * Film profile loading from JSON
 *
 * Uses GLib's json-glib parser
 */

#include "filmsim_common.h"
#include <string.h>
#include <stdio.h>
#include <glib.h>
#include <json-glib/json-glib.h>

bool filmsim_load_profile_from_json(filmsim_film_profile_t *profile,
                                     const char *json_string) {
    JsonParser *parser = NULL;
    JsonNode *root = NULL;
    JsonObject *obj = NULL;
    bool success = false;

    /* Initialize profile to zero */
    memset(profile, 0, sizeof(filmsim_film_profile_t));
    
    /* Set default dye mix matrix to identity (no mixing) */
    profile->couplers.dye_mix_matrix[0][0] = 1.0;
    profile->couplers.dye_mix_matrix[1][1] = 1.0;
    profile->couplers.dye_mix_matrix[2][2] = 1.0;

    /* Parse JSON */
    parser = json_parser_new();
    GError *error = NULL;

    if (!json_parser_load_from_data(parser, json_string, -1, &error)) {
        g_error_free(error);
        goto cleanup;
    }

    root = json_parser_get_root(parser);
    if (!root || !JSON_NODE_HOLDS_OBJECT(root)) {
        goto cleanup;
    }

    obj = json_node_get_object(root);

    /* Parse info section */
    if (json_object_has_member(obj, "info")) {
        JsonObject *info_obj = json_object_get_object_member(obj, "info");

        const char *name = json_object_get_string_member(info_obj, "name");
        if (name) {
            strncpy(profile->info.name, name, sizeof(profile->info.name) - 1);
            profile->info.name[sizeof(profile->info.name) - 1] = '\0';
        }

        const char *desc = json_object_get_string_member(info_obj, "description");
        if (desc) {
            strncpy(profile->info.description, desc, sizeof(profile->info.description) - 1);
            profile->info.description[sizeof(profile->info.description) - 1] = '\0';
        }

        if (json_object_has_member(info_obj, "format_mm")) {
            profile->info.format_mm = json_object_get_int_member(info_obj, "format_mm");
        }
        
        /* Parse film_type */
        if (json_object_has_member(info_obj, "film_type")) {
            const char *film_type = json_object_get_string_member(info_obj, "film_type");
            if (film_type) {
                if (g_strcmp0(film_type, "monochrome_negative") == 0) {
                    profile->info.film_type = FILMSIM_FILM_TYPE_MONOCHROME_NEGATIVE;
                } else if (g_strcmp0(film_type, "color_positive_reversal") == 0 ||
                           g_strcmp0(film_type, "color_positive") == 0) {
                    profile->info.film_type = FILMSIM_FILM_TYPE_COLOR_POSITIVE;
                } else {
                    profile->info.film_type = FILMSIM_FILM_TYPE_COLOR_NEGATIVE;
                }
            }
        }
    }

    /* Parse processing section (unused in plugin but load for completeness) */
    if (json_object_has_member(obj, "processing")) {
        (void)json_object_get_object_member(obj, "processing"); /* unused */
        /* Just store, not used in pipeline */
    }

    /* Parse properties section */
    if (!json_object_has_member(obj, "properties")) {
        goto cleanup;
    }
    JsonObject *props = json_object_get_object_member(obj, "properties");

    /* Parse calibration */
    if (json_object_has_member(props, "calibration")) {
        JsonObject *cal_obj = json_object_get_object_member(props, "calibration");
        profile->calibration.iso = json_object_get_int_member(cal_obj, "iso");
        profile->calibration.middle_gray_logE = json_object_get_double_member(cal_obj, "middle_gray_logE");
    }

    /* Parse halation */
    if (json_object_has_member(props, "halation")) {
        JsonObject *hal_obj = json_object_get_object_member(props, "halation");

        if (json_object_has_member(hal_obj, "strength")) {
            JsonObject *str_obj = json_object_get_object_member(hal_obj, "strength");
            profile->halation.strength.r = json_object_get_double_member(str_obj, "r");
            profile->halation.strength.g = json_object_get_double_member(str_obj, "g");
            profile->halation.strength.b = json_object_get_double_member(str_obj, "b");
        }

        if (json_object_has_member(hal_obj, "size_um")) {
            JsonObject *size_obj = json_object_get_object_member(hal_obj, "size_um");
            profile->halation.size_um.r = json_object_get_double_member(size_obj, "r");
            profile->halation.size_um.g = json_object_get_double_member(size_obj, "g");
            profile->halation.size_um.b = json_object_get_double_member(size_obj, "b");
        }
    }

    /* Parse couplers */
    if (json_object_has_member(props, "couplers")) {
        JsonObject *couplers_obj = json_object_get_object_member(props, "couplers");

        profile->couplers.saturation_amount = json_object_get_double_member(couplers_obj, "saturation_amount");

        if (json_object_has_member(couplers_obj, "dir_amount_rgb")) {
            JsonArray *dir_array = json_object_get_array_member(couplers_obj, "dir_amount_rgb");
            guint len = json_array_get_length(dir_array);
            for (guint i = 0; i < len && i < 3; i++) {
                profile->couplers.dir_amount_rgb[i] = json_array_get_double_element(dir_array, i);
            }
        }

        profile->couplers.dir_diffusion_um = json_object_get_double_member(couplers_obj, "dir_diffusion_um");
        profile->couplers.dir_diffusion_interlayer = json_object_get_double_member(couplers_obj, "dir_diffusion_interlayer");

        /* Parse dye mixing matrix */
        if (json_object_has_member(couplers_obj, "dye_mix_matrix")) {
            JsonArray *matrix_array = json_object_get_array_member(couplers_obj, "dye_mix_matrix");
            guint rows = json_array_get_length(matrix_array);
            for (guint i = 0; i < rows && i < 3; i++) {
                JsonArray *row_array = json_array_get_array_element(matrix_array, i);
                guint cols = json_array_get_length(row_array);
                for (guint j = 0; j < cols && j < 3; j++) {
                    profile->couplers.dye_mix_matrix[i][j] = json_array_get_double_element(row_array, j);
                }
            }
        }
    }

    /* Parse curves (H-D curves) */
    if (json_object_has_member(props, "curves")) {
        JsonObject *curves_obj = json_object_get_object_member(props, "curves");

        if (json_object_has_member(curves_obj, "hd")) {
            JsonArray *hd_array = json_object_get_array_member(curves_obj, "hd");
            guint len = json_array_get_length(hd_array);
            profile->num_curves = (len < FILMSIM_MAX_CURVE_POINTS) ? len : FILMSIM_MAX_CURVE_POINTS;

            for (guint i = 0; i < profile->num_curves; i++) {
                JsonNode *pt_node = json_array_get_element(hd_array, i);
                if (!JSON_NODE_HOLDS_OBJECT(pt_node)) continue;

                JsonObject *pt_obj = json_node_get_object(pt_node);
                profile->curves.rgb[i].d = json_object_get_double_member(pt_obj, "d");
                
                /* Check if this is RGB curves or mono (density) curves */
                if (json_object_has_member(pt_obj, "r")) {
                    /* RGB curves (color film) */
                    profile->curves.rgb[i].r = json_object_get_double_member(pt_obj, "r");
                    profile->curves.rgb[i].g = json_object_get_double_member(pt_obj, "g");
                    profile->curves.rgb[i].b = json_object_get_double_member(pt_obj, "b");
                } else if (json_object_has_member(pt_obj, "density")) {
                    /* Mono curves (B&W film) - use same density for all channels */
                    double density = json_object_get_double_member(pt_obj, "density");
                    profile->curves.rgb[i].r = density;
                    profile->curves.rgb[i].g = density;
                    profile->curves.rgb[i].b = density;
                }
            }
        }

        /* Parse spectral sensitivity (color films: y, m, c; mono: sensitivity) */
        if (json_object_has_member(curves_obj, "spectral_sensitivity")) {
            JsonArray *ss_array = json_object_get_array_member(curves_obj, "spectral_sensitivity");
            guint len = json_array_get_length(ss_array);
            profile->num_spectral_sensitivity = (len < FILMSIM_MAX_SPECTRAL_BANDS) ? len : FILMSIM_MAX_SPECTRAL_BANDS;

            for (guint i = 0; i < profile->num_spectral_sensitivity; i++) {
                JsonNode *pt_node = json_array_get_element(ss_array, i);
                if (!JSON_NODE_HOLDS_OBJECT(pt_node)) continue;

                JsonObject *pt_obj = json_node_get_object(pt_node);
                double wavelength = json_object_get_double_member(pt_obj, "wavelength");

                if (json_object_has_member(pt_obj, "y")) {
                    /* Color film: Y, M, C layer sensitivities */
                    profile->spectral_sensitivity.rgb[i].wavelength = wavelength;
                    profile->spectral_sensitivity.rgb[i].y = json_object_get_double_member(pt_obj, "y");
                    profile->spectral_sensitivity.rgb[i].m = json_object_get_double_member(pt_obj, "m");
                    profile->spectral_sensitivity.rgb[i].c = json_object_get_double_member(pt_obj, "c");
                } else if (json_object_has_member(pt_obj, "sensitivity")) {
                    /* Mono film: single sensitivity value */
                    profile->spectral_sensitivity.mono[i].wavelength = wavelength;
                    profile->spectral_sensitivity.mono[i].sensitivity = json_object_get_double_member(pt_obj, "sensitivity");
                }
            }
        }

        /* Parse spectral density at D-min (film base, for negative films) */
        if (json_object_has_member(curves_obj, "spectral_density_min")) {
            JsonArray *sd_array = json_object_get_array_member(curves_obj, "spectral_density_min");
            guint len = json_array_get_length(sd_array);
            profile->num_spectral_density_min = (len < FILMSIM_MAX_SPECTRAL_BANDS) ? len : FILMSIM_MAX_SPECTRAL_BANDS;

            for (guint i = 0; i < profile->num_spectral_density_min; i++) {
                JsonNode *pt_node = json_array_get_element(sd_array, i);
                if (!JSON_NODE_HOLDS_OBJECT(pt_node)) continue;

                JsonObject *pt_obj = json_node_get_object(pt_node);
                profile->spectral_density_min[i].wavelength = json_object_get_double_member(pt_obj, "wavelength");
                profile->spectral_density_min[i].density = json_object_get_double_member(pt_obj, "density");
            }
        }

        /* Parse spectral density at D-mid (mid-gray exposure, for negative films) */
        if (json_object_has_member(curves_obj, "spectral_density_mid")) {
            JsonArray *sd_array = json_object_get_array_member(curves_obj, "spectral_density_mid");
            guint len = json_array_get_length(sd_array);
            profile->num_spectral_density_mid = (len < FILMSIM_MAX_SPECTRAL_BANDS) ? len : FILMSIM_MAX_SPECTRAL_BANDS;

            for (guint i = 0; i < profile->num_spectral_density_mid; i++) {
                JsonNode *pt_node = json_array_get_element(sd_array, i);
                if (!JSON_NODE_HOLDS_OBJECT(pt_node)) continue;

                JsonObject *pt_obj = json_node_get_object(pt_node);
                profile->spectral_density_mid[i].wavelength = json_object_get_double_member(pt_obj, "wavelength");
                profile->spectral_density_mid[i].density = json_object_get_double_member(pt_obj, "density");
            }
        }

        /* Parse per-dye spectral density curves (for reversal/slide films) */
        /* Yellow dye */
        if (json_object_has_member(curves_obj, "spectral_density_y")) {
            JsonArray *sd_array = json_object_get_array_member(curves_obj, "spectral_density_y");
            guint len = json_array_get_length(sd_array);
            profile->num_spectral_density_y = (len < FILMSIM_MAX_SPECTRAL_BANDS) ? len : FILMSIM_MAX_SPECTRAL_BANDS;

            for (guint i = 0; i < profile->num_spectral_density_y; i++) {
                JsonNode *pt_node = json_array_get_element(sd_array, i);
                if (!JSON_NODE_HOLDS_OBJECT(pt_node)) continue;

                JsonObject *pt_obj = json_node_get_object(pt_node);
                profile->spectral_density_y[i].wavelength = json_object_get_double_member(pt_obj, "wavelength");
                profile->spectral_density_y[i].density = json_object_get_double_member(pt_obj, "density");
            }
        }

        /* Magenta dye */
        if (json_object_has_member(curves_obj, "spectral_density_m")) {
            JsonArray *sd_array = json_object_get_array_member(curves_obj, "spectral_density_m");
            guint len = json_array_get_length(sd_array);
            profile->num_spectral_density_m = (len < FILMSIM_MAX_SPECTRAL_BANDS) ? len : FILMSIM_MAX_SPECTRAL_BANDS;

            for (guint i = 0; i < profile->num_spectral_density_m; i++) {
                JsonNode *pt_node = json_array_get_element(sd_array, i);
                if (!JSON_NODE_HOLDS_OBJECT(pt_node)) continue;

                JsonObject *pt_obj = json_node_get_object(pt_node);
                profile->spectral_density_m[i].wavelength = json_object_get_double_member(pt_obj, "wavelength");
                profile->spectral_density_m[i].density = json_object_get_double_member(pt_obj, "density");
            }
        }

        /* Cyan dye */
        if (json_object_has_member(curves_obj, "spectral_density_c")) {
            JsonArray *sd_array = json_object_get_array_member(curves_obj, "spectral_density_c");
            guint len = json_array_get_length(sd_array);
            profile->num_spectral_density_c = (len < FILMSIM_MAX_SPECTRAL_BANDS) ? len : FILMSIM_MAX_SPECTRAL_BANDS;

            for (guint i = 0; i < profile->num_spectral_density_c; i++) {
                JsonNode *pt_node = json_array_get_element(sd_array, i);
                if (!JSON_NODE_HOLDS_OBJECT(pt_node)) continue;

                JsonObject *pt_obj = json_node_get_object(pt_node);
                profile->spectral_density_c[i].wavelength = json_object_get_double_member(pt_obj, "wavelength");
                profile->spectral_density_c[i].density = json_object_get_double_member(pt_obj, "density");
            }
        }
    }

    success = true;

cleanup:
    if (parser) g_object_unref(parser);
    return success;
}

/*
 * Static storage for loaded profile JSON
 * (Avoids memory management complexity - single profile at a time)
 */
static gchar *loaded_profile_json = NULL;

/*
 * Search paths for film profiles
 */
static const char *profile_search_paths[] = {
    "../data",                                  /* Relative to plugin */
    "../../data",                               /* When in build dir */
    "/usr/share/cppfilmsim/data",               /* System-wide install */
    "/usr/local/share/cppfilmsim/data",         /* Local install */
    NULL
};

/*
 * Get embedded film profile by name
 * Returns JSON string or NULL if not found
 * The returned string is valid until next call to this function
 */
const char* filmsim_get_embedded_profile(const char *film_name) {
    if (!film_name) return NULL;
    
    /* Free previous profile */
    if (loaded_profile_json) {
        g_free(loaded_profile_json);
        loaded_profile_json = NULL;
    }
    
    char filename[256];
    snprintf(filename, sizeof(filename), "%s.json", film_name);
    
    /* Try each search path */
    for (int i = 0; profile_search_paths[i] != NULL; i++) {
        char path[512];
        snprintf(path, sizeof(path), "%s/%s", profile_search_paths[i], filename);
        
        gsize length = 0;
        if (g_file_get_contents(path, &loaded_profile_json, &length, NULL)) {
            return loaded_profile_json;
        }
    }
    
    /* Try user config directory */
    const char *config_dir = g_get_user_config_dir();
    if (config_dir) {
        char path[512];
        snprintf(path, sizeof(path), "%s/darktable/filmsim/%s", config_dir, filename);
        
        gsize length = 0;
        if (g_file_get_contents(path, &loaded_profile_json, &length, NULL)) {
            return loaded_profile_json;
        }
    }
    
    return NULL;
}

/*
 * List all embedded film profiles
 * Returns number of profiles, fills names array with allocated strings
 * Caller must free the strings and array
 */
int filmsim_list_embedded_profiles(char ***names) {
    static const char *film_names[] = {
        "kodak_portra_400",
        "kodak_ektar_100",
        "kodak_ultramax_400",
        "kodak_vision3_500t",
        "kodak_kodachrome_64",
        "kodak_ektachrome_p_1600",
        "kodak_royal_gold_1000",
        "kodak_portra_100T",
        "kodak_t-max_400_tmax-rs_7min",
        "fuji_proh_400",
        "fuji_superia_1600",
        "fuji_superia_reala_100",
        "fuji_fujichrome_velvia_100",
        "fuji_fujichrome_astia_f_100",
        "fuji_neopan_arcos_ii_100",
        "cinestill_800t",
        "agfacolor_optima_ii_100",
        "agfacolor_portrait_xps_160",
        NULL
    };
    
    int count = 0;
    while (film_names[count]) count++;
    
    if (names) {
        *names = (char **)g_malloc((count + 1) * sizeof(char *));
        for (int i = 0; i < count; i++) {
            (*names)[i] = g_strdup(film_names[i]);
        }
        (*names)[count] = NULL;
    }
    
    return count;
}
