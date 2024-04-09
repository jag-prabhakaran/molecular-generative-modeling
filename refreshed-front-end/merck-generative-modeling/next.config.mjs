/** @type {import('next').NextConfig} */

import CopyPlugin from "copy-webpack-plugin";

const nextConfig = {
  output: "export",
  webpack: (config) => {
    config.plugins.push(
      new CopyPlugin({
        patterns: [{ from: "public", to: "" }],
      })
    );
    return config;
  },
};

export default nextConfig;
