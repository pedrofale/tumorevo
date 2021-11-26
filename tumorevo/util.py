import numpy as np
import pandas as pd

def prepare_muller(genotype_counts, genotype_parents):
        pop_df = pd.DataFrame(genotype_counts).fillna(0)
        pop_df['Generation'] = np.arange(pop_df.shape[0])
        pop_df = pop_df.melt(id_vars=['Generation'], value_name='Population', var_name='Identity').sort_values('Generation').reset_index(drop=True)

        anc_df = pd.DataFrame([genotype_parents[-1]]).melt(var_name='Identity', value_name='Parent')

        color_by = pd.Series(np.arange(anc_df.shape[0]+1), index=pop_df['Identity'].unique())

        return pop_df, anc_df, color_by
