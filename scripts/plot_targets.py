from dna_features_viewer import (
    GraphicFeature,
    GraphicRecord,
    CircularGraphicRecord,
    BiopythonTranslator,
)
from colour import Color
import yaml

Color('blue')

def read_relabel(filepath):
    with open(filepath) as handle:
        try:
            return yaml.safe_load(handle)
        except yaml.YAMLError as exc:
            print(exc)


def relabel(features, relabel_dict):
    """Remove / edit feature labels so resulting plot looks cleaner.

    Args:
        features ([list]): List of features.
    """
    delete_indexes = []
    for each_relab in relabel_dict:
        feat_relab_dict = relabel_dict[each_relab]
        for i, each_feat in enumerate(features):
            if feat_relab_dict["search"] in each_feat.label:  # search term in label
                if "remove" in feat_relab_dict:  # if relabel dict has remove key
                    delete_indexes.append(i)     # mark feature for destruction
                else:  # otherwise update feature attributes using values in
                       # dictionary stored under attrs
                    features[i].__dict__.update(feat_relab_dict["attrs"])
                    if features[i].label and features[i].label == "":
                        features[i].label = "test"

    drawn_features = [feat for i, feat in enumerate(features) if i not in delete_indexes]
    return drawn_features


def relabel_cas9_targets(features):
    """Relabel and recolor Cas9 target sites with nice gradient.

    Args:
        features (list): Feature list.
    """
    # get Cas9 target site features indexes
    targets = [i for i, feat in enumerate(features) if "Cas9" in feat.label]
    start_color = Color("#73FBD3")
    end_color = Color("#5C7AFF")
    gradient = [c.hex for c in list(start_color.range_to(end_color, steps=len(targets)))]
    for i, each_index in enumerate(targets):
        features[each_index].label = f"Cas9 T{i}"
        features[each_index].color = gradient[i]
    return features


def draw_map(record_path, relabel_dict, output_path):

    gr = BiopythonTranslator().translate_record(record_path)
    feats = gr.features
    feats = relabel(feats, relabel_dict)
    feats = relabel_cas9_targets(feats)
    #print(len(feats))
    circ_gr = CircularGraphicRecord(sequence_length=len(gr.sequence), features=feats)
    ax, _ = circ_gr.plot()
    ax.figure.savefig(output_path, bbox_inches='tight')


def main():

    relabel_dict = read_relabel(snakemake.input["relabel_dict"])['relabel']
    draw_map(str(snakemake.input["record"]), relabel_dict, str(snakemake.output))


if __name__ == "__main__":
    main()
