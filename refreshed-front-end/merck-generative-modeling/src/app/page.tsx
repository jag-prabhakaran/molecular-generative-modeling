import MerckNavbar from "@/app/_components/MerckNavbar";
import ModelCard from "@/app/_components/ModelCard";
import { Box, CssBaseline, Toolbar } from "@mui/material";

export default function Home() {
  const models = [];
  for (let i = 0; i < 7; i++) {
    models.push(<ModelCard />);
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
