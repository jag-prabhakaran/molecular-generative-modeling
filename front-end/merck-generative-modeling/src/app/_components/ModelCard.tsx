import * as React from "react";
import Box from "@mui/material/Box";
import Card from "@mui/material/Card";
import CardActions from "@mui/material/CardActions";
import CardContent from "@mui/material/CardContent";
import Button from "@mui/material/Button";
import Typography from "@mui/material/Typography";

const bull = (
  <Box
    component="span"
    sx={{ display: "inline-block", mx: "2px", transform: "scale(0.8)" }}
  >
    •
  </Box>
);

export default function ModelCard(props: {
  ModelTitle: string;
  ModelDescription: string;
  ModelLink: string;
}) {
  // Now you can destructure with types already defined
  const { ModelTitle, ModelDescription, ModelLink } = props;

  return (
    <Box sx={{ minWidth: 275 }}>
      <Card variant="outlined">
        <React.Fragment>
          <CardContent>
            <Typography
              sx={{ fontSize: 14 }}
              color="text.secondary"
              gutterBottom
            >
              {ModelTitle} {/* Use the ModelTitle prop */}
            </Typography>
            <Typography variant="body2">
              {ModelDescription} {/* Use the ModelDescription prop */}
            </Typography>
          </CardContent>
          <CardActions>
            <Button size="small" href={ModelLink}>
              Run Model
            </Button>
          </CardActions>
        </React.Fragment>
      </Card>
    </Box>
  );
}
