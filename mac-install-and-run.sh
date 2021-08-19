echo "Installing Homebrew"
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

echo "Installing CMake"
brew install cmake

echo "Installing Eigen"
brew install eigen

echo "Building StrokeStrip"
mkdir build
cd build
cmake ..
make
cd ..

echo "Generating content for Fig. 17"
./examples/letters/generate_letters_data.sh

echo "Done!"
echo "See the outputted SVG files in examples/letters"
