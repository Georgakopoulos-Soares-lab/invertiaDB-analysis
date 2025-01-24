from setuptools import setup, find_packages


requirements = []
with open("requirements.txt") as f:
    for line in f:
        line = line.strip()
        requirements.append(line)

setup(
    name='mindi',
    version='0.1',                
    packages=find_packages(),    
    description='MINDI: Your favourite non B-DNA package',
    author='Nikol Chantzi',
    author_email='nicolechantzi@gmail.com',
    url='https://github.com/yourusername/my_package',  
    install_requires=requirements
)
