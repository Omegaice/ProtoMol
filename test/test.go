package main

import (
	"bufio"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
)

func main() {
	log.SetFlags(log.LstdFlags | log.Lshortfile)

	// Update Path Variable
	os.Setenv("PATH", ".:"+os.Getenv("PATH"))

	// Check if we can find ProtoMol
	if _, err := exec.LookPath("ProtoMol"); err != nil {
		log.Fatal("ProtoMol executable could not be found")
	}

	// Find Tests
	tests, err := filepath.Glob("tests/*.conf")
	if err != nil {
		log.Fatalln(err)
	}

	// Create Output Directory
	os.Mkdir("tests/output", 0755)

	// Run Tests
	for _, test := range tests {
		basename := strings.Replace(filepath.Base(test), filepath.Ext(test), "", -1)
		log.Println("Executing Test:", basename)

		cmd := exec.Command("ProtoMol", test)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			log.Fatal(err)
		}

		// Find Expected Results
		expects, err := filepath.Glob("tests/expected/" + basename + ".*")
		if err != nil {
			log.Fatalln(err)
		}

		// Compare Results
		for _, expected := range expects {
			output := "tests/output/" + filepath.Base(expected)

			extension := filepath.Ext(output)
			log.Println(extension)
			switch extension {
			case ".dcd":
				if !isMatchingDCD(output, expected) {
					log.Println("Differ")
				}
				break
			case ".energy":
				if !isMatchingEnergy(output, expected) {
					log.Println("Differ")
				}
				break
			case ".forces":
				break
			case ".pos":
				break
			case ".vel":
				break
			}
		}

		break
	}
}

// DCD Comparison
func isMatchingDCD(actual, expected string) bool {
	dcdActual, err := ReadDCD(actual)
	if err != nil {
		return true
	}

	dcdExpected, err := ReadDCD(expected)
	if err != nil {
		return true
	}

	if dcdActual.Header.Frames != dcdExpected.Header.Frames {
		log.Println("Atom Count Differs")
		return false
	}

	if dcdActual.Header.FirstStep != dcdExpected.Header.FirstStep {
		log.Printf("First Step Differs. Should be %d but is %d\n", dcdExpected.Header.FirstStep, dcdActual.Header.FirstStep)
		return false
	}

	if dcdActual.Atoms != dcdExpected.Atoms {
		log.Printf("Atom Count Differs. Should be %d but is %d\n", dcdExpected.Atoms, dcdActual.Atoms)
		return false
	}

	diffs := 0
	for frame := int32(0); frame < dcdExpected.Header.Frames; frame++ {
		for atom := int32(0); atom < dcdExpected.Atoms; atom++ {
			xExpected := dcdExpected.Frames[frame].X[atom]
			xActual := dcdActual.Frames[frame].X[atom]

			if math.Max(float64(xExpected), float64(xActual))-math.Min(float64(xExpected), float64(xActual)) > 0.00001 {
				diffs++
				log.Printf("Frame %d, Atom %d Differs\n", frame, atom)
				log.Printf("Expected: %f, Actual: %f, Difference: %f\n", xExpected, xActual, math.Max(float64(xExpected), float64(xActual))-math.Min(float64(xExpected), float64(xActual)))
			}

			yExpected := dcdExpected.Frames[frame].Y[atom]
			yActual := dcdActual.Frames[frame].Y[atom]

			if math.Max(float64(yExpected), float64(yActual))-math.Min(float64(yExpected), float64(yActual)) > 0.00001 {
				diffs++
				log.Printf("Frame %d, Atom %d Differs\n", frame, atom)
				log.Printf("Expected: %f, Actual: %f, Difference: %f\n", yExpected, yActual, math.Max(float64(yExpected), float64(yActual))-math.Min(float64(yExpected), float64(yActual)))
			}

			zExpected := dcdExpected.Frames[frame].X[atom]
			zActual := dcdActual.Frames[frame].Z[atom]

			if math.Max(float64(zExpected), float64(zActual))-math.Min(float64(zExpected), float64(zActual)) > 0.00001 {
				diffs++
				log.Printf("Frame %d, Atom %d Differs\n", frame, atom)
				log.Printf("Expected: %f, Actual: %f, Difference: %f\n", zExpected, zActual, math.Max(float64(zExpected), float64(zActual))-math.Min(float64(zExpected), float64(zActual)))
			}
		}
	}

	return diffs == 0
}

// DCD Parsing
type DCDFile struct {
	Header  DCDHeader
	Comment []byte
	Atoms   int32
	Frames  []DCDFrame
}

type DCDFrame struct {
	X, Y, Z []float32
}
type DCDHeader struct {
	Frames      int32
	FirstStep   int32
	_           [24]byte
	FreeIndexes int32
	_           [44]byte
	_           [4]byte
}

func ReadDCD(filename string) (*DCDFile, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	// Read Endian
	var endian int32
	if err := binary.Read(file, binary.LittleEndian, &endian); err != nil {
		return nil, err
	}

	// Test Endian
	var order binary.ByteOrder = binary.LittleEndian
	if endian != 84 {
		order = binary.BigEndian
	}

	// Read Cord
	var cord [4]byte
	if err := binary.Read(file, order, &cord); err != nil {
		return nil, err
	}

	// Check Magic Number
	if cord != [4]byte{'C', 'O', 'R', 'D'} {
		return nil, err
	}

	var result DCDFile

	// Read Header
	if err := binary.Read(file, order, &result.Header); err != nil {
		return nil, err
	}

	// Read Comment
	if err := FortranReadString(file, order, &result.Comment); err != nil {
		return nil, err
	}

	// Read Atoms
	if err := FortranRead(file, order, &result.Atoms); err != nil {
		return nil, err
	}

	// Skip Free Indexies
	if result.Header.FreeIndexes > 0 {
		io.CopyN(ioutil.Discard, file, int64(4*(result.Atoms-result.Header.FreeIndexes+2)))
	}

	// Read Frames
	for i := int32(0); i < result.Header.Frames; i++ {
		var frame DCDFrame

		// Read X
		frame.X = make([]float32, result.Atoms)
		if err := FortranRead(file, order, &frame.X); err != nil {
			return nil, err
		}

		// Read Y
		frame.Y = make([]float32, result.Atoms)
		if err := FortranRead(file, order, &frame.Y); err != nil {
			return nil, err
		}

		// Read Z
		frame.Z = make([]float32, result.Atoms)
		if err := FortranRead(file, order, &frame.Z); err != nil {
			return nil, err
		}

		result.Frames = append(result.Frames, frame)
	}

	return &result, nil
}

func FortranRead(file io.Reader, order binary.ByteOrder, data interface{}) error {
	var head int32
	if err := binary.Read(file, order, &head); err != nil {
		return err
	}

	if err := binary.Read(file, order, data); err != nil {
		return err
	}

	var tail int32
	if err := binary.Read(file, order, &tail); err != nil {
		return err
	}

	if head != tail {
		return errors.New(fmt.Sprintf("head (%d) does not match tail (%d)", head, tail))
	}

	return nil
}

func FortranReadString(file io.Reader, order binary.ByteOrder, data *[]byte) error {
	var head int32
	if err := binary.Read(file, order, &head); err != nil {
		return err
	}

	*data = make([]byte, head)
	if err := binary.Read(file, order, data); err != nil {
		return err
	}

	var tail int32
	if err := binary.Read(file, order, &tail); err != nil {
		return err
	}

	return nil
}

// Energy Comparison
func isMatchingEnergy(actual, expected string) bool {
	energyActual, err := ReadEnergy(actual)
	if err != nil {
		return true
	}

	energyExpected, err := ReadEnergy(expected)
	if err != nil {
		return true
	}

	if len(energyActual.Column) != len(energyExpected.Column) {
		log.Printf("Column Count Differs. Should be %d but is %d\n", len(energyActual.Column), len(energyExpected.Column))
		return false
	}

	var diffs int
	for i := 0; i < len(energyExpected.Column); i++ {
		if energyActual.Column[i].Name != energyExpected.Column[i].Name {
			log.Printf("Column Name Differs. Should be %s but is %s\n", energyActual.Column[i].Name, energyExpected.Column[i].Name)
			return false
		}

		if len(energyActual.Column[i].Values) != len(energyExpected.Column[i].Values) {
			log.Printf("Column Value Count Differs. Should be %d but is %d\n", len(energyActual.Column[i].Values), len(energyExpected.Column[i].Values))
			return false
		}

		for j := 0; j < len(energyExpected.Column[i].Values); j++ {
			xExpected := energyExpected.Column[i].Values[j]
			xActual := energyActual.Column[i].Values[j]

			if math.Max(float64(xExpected), float64(xActual))-math.Min(float64(xExpected), float64(xActual)) > 0.00001 {
				diffs++
				log.Printf("Column Value Differs. Expected: %f, Actual: %f, Difference: %f\n", xExpected, xActual, math.Max(float64(xExpected), float64(xActual))-math.Min(float64(xExpected), float64(xActual)))
			}
		}
	}

	return diffs == 0
}

// Energy Parsing
type EnergyFile struct {
	Column []EnergyFileColumn
}

type EnergyFileColumn struct {
	Name   string
	Values []float32
}

func ReadEnergy(filename string) (*EnergyFile, error) {
	// Parse Headers
	sHeader, err := ioutil.ReadFile(fmt.Sprintf("%s.header", filename))
	if err != nil {
		return nil, err
	}

	var result EnergyFile

	for _, header := range strings.Split(string(sHeader), " ") {
		if len(header) == 0 {
			continue
		}
		result.Column = append(result.Column, EnergyFileColumn{Name: header})
	}

	// Parse File
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		columnID := 0
		for _, column := range strings.Split(scanner.Text(), " ") {
			if len(column) == 0 {
				continue
			}

			value, err := strconv.ParseFloat(column, 32)
			if err != nil {
				return nil, err
			}

			result.Column[columnID].Values = append(result.Column[columnID].Values, float32(value))
			columnID++
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return &result, nil
}
