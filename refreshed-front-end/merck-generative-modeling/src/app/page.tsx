import MerckNavbar from "@/app/_components/MerckNavbar";
import ModelCard from "@/app/_components/ModelCard";
import { Mode } from "@mui/icons-material";
import { Box, CssBaseline, Toolbar } from "@mui/material";

export default function Home() {
  return (
    <Box className="flex">
      <CssBaseline />
      <MerckNavbar />
      <Box component="main" className="flex flex-col flex-grow">
        <Toolbar />
        <ModelCard
          ModelTitle="VAE GAN"
          ModelDescription="This model is a simple GAN (Generative Adversarial Network) architecture that trains on small molecules (number of atoms < 9) and outputs new (small) molecules from noise."
          ModelLink="./models/vae-gan"
        />
        <ModelCard
          ModelTitle="SAMOA"
          ModelDescription="This model employs a SMILES-based Recurrent Neural Network (RNN) generative model to 
        achieve scaffold-constrained generation."
          ModelLink="./models/scaffold-constrained"
        />
        <ModelCard
          ModelTitle="Graph Translation"
          ModelDescription="The model works by decomposing a molecule into a graph where each node corresponds to a specific substructure 
        of the molecule. This graph is then passed through a VAE (Variational Autoencoder) model to generate new molecules."
          ModelLink="./models/multiobj-rationale"
        />
      </Box>
    </Box>
  );
}
