import MerckNavbar from "@/app/_components/MerckNavbar";
import ModelCard from "@/app/_components/ModelCard";
import { Box, CssBaseline, Toolbar } from "@mui/material";

export default function Home() {
  const models = [];
  for (let i = 0; i < 3; i++) {
    models.push(<ModelCard 
      ModelTitle="VAE GAN"
      ModelDescription="simple GAN (Generative Adversarial Network) architecture that trains on small molecules (number of atoms < 9) and outputs new (small) molecules from noise."
    />);
  }
  return (
    <Box className="flex">
      <CssBaseline />
      <MerckNavbar />
      <Box component="main" className="flex flex-col flex-grow">
        <Toolbar />
        {models}
      </Box>
    </Box>
  );
}
